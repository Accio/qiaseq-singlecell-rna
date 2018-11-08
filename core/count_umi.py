import itertools
import logging
import sys
from guppy import hpy
import pysam
from pathos.multiprocessing import ProcessingPool as Pool
import regex
from intervaltree import IntervalTree
from collections import defaultdict,OrderedDict
from functools import partial
from Bio.Seq import Seq

## Modules from this project
from find_primer import find_primer
from find_gene import find_gene
from demultiplex_cells import write_metrics
from create_annotation_tables import create_gene_tree

def grouper(iterable,n=750000):
    '''
    Returns n chunks of an iterable

    :param iterator iterable: The iterator to chop up
    :param int n: the chunks to group the iterator into
    :yields: iterator , i.e. the n1, n2, ... n chunks of the iterable
    '''
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it,n))
        if not chunk:
            return
        yield chunk

def iterate_bam_chunks(tagged_bam,chunks=750000):
    '''
    Iterate over a bam file in chunks (i.e. the number of reads returned
    at a time)

    :param str tagged_bam: the input bam file with UMI tags
    :param int chunks: the number of reads to process at a time
    :yields: a list of tuples
    '''
    IN = pysam.AlignmentFile(tagged_bam,'rb')
    chroms = IN.header['SQ']
    for reads in grouper(IN.fetch(until_eof=True),chunks):
        to_yield = []
        for read in reads:
            if read.flag == 4: ## Unmapped
                chromosome = '*'
            else:
                chromosome = chroms[read.tid]['SN']
            umi = read.qname.split(":")[-1]
            to_yield.append((read.qname,read.seq, read.is_reverse, read.alen,chromosome,
                             read.pos, read.cigarstring, umi,
                             read.get_tag('NH')))
        yield to_yield

def count_umis_wts(gene_tree,tagged_bam,outfile,metricfile,logfile,cores=3):
    ''' Count UMIs for each gene in the input the tagged_bam file

    :param object gene_tree : an IntervalTree data structure
    :param str tagged_bam: a UMI tagged bam file
    :param str outfile: the output file
    :param str metricfile: file to write the metrics stats
    :param str logfile: the log file to write to
    :param int cores: the number of cores to use
    '''
    ## Set up logging
    logger = logging.getLogger("count_umis")
    logger.setLevel(logging.DEBUG)
    LOG = logging.FileHandler(logfile)
    logger.addHandler(LOG)
    ## Variable Initialization
    found = 0
    found_ercc = 0
    miss_chr = 0
    unmapped = 0
    not_annotated = 0
    multimapped = 0
    multimapped_ercc = 0
    ercc  = 0
    total_UMIs = 0
    ## Store umi counts for each gene
    umi_counter = defaultdict(lambda:defaultdict(int))
    logger.info('Using {} cores'.format(cores))
    p = Pool(cores)
    func = partial(find_gene,gene_tree)
    max_reads_in_mem = 10000000
    
    for reads in iterate_bam_chunks(tagged_bam,chunks=max_reads_in_mem):
        logger.info('Reading {} reads in memory to find genes'.format(max_reads_in_mem))
        find_gene_results = p.map(func,reads)
        for info in find_gene_results:
            gene_info,umi,count,nh = info
            if count == 0:
                if gene_info == 'Unknown_Chrom':
                    miss_chr+=1
                elif gene_info == 'Unmapped':
                    unmapped+=1
                elif gene_info == 'Unknown':
                    not_annotated+=1
                else: ## These are ERCC reads
                    ercc+=1
                    if nh>1:
                        multimapped_ercc+=1
                        multimapped+=1
                    else:
                        umi_counter[gene_info][umi]+=1
                        found+=1
                        found_ercc+=1
            else:
                if nh>1:
                    multimapped+=1
                else:
                    umi_counter[gene_info][umi]+=1
                    found+=1
    p.close()
    p.join()
    ## Print output results
    ## Write gene counts
    detected_genes = set()
    with open(outfile,'w') as OUT:
        for chrom in gene_tree:
            for strand in gene_tree[chrom]:
                for gene_info in gene_tree[chrom][strand]:
                    ensembl_id,gene,strand,chrom,five_prime,three_prime = gene_info.data
                    if gene_info.data in umi_counter:
                        umi_count = len(umi_counter[gene_info.data])
                        if not gene.startswith('ERCC'):
                            detected_genes.add(gene_info.data)
                    else:
                        umi_count = 0
                    total_UMIs+=umi_count
                    OUT.write(ensembl_id+'\t'+gene+"\t"+strand+"\t"+chrom+"\t"+str(five_prime)+'\t'+str(three_prime)+"\t"+str(umi_count)+"\n")
    ## Write metrics
    metric_dict = OrderedDict([
        ('reads dropped, not mapped to genome',unmapped),
        ('reads dropped, not annotated',not_annotated+miss_chr),
        ('reads dropped, aligned to genome, multiple loci',multimapped-multimapped_ercc),
        ('reads dropped, aligned to ERCC, multiple loci',multimapped_ercc),        
        ('reads used, aligned to genome, unique loci',found-found_ercc),
        ('reads used, aligned to ERCC, unique loci',found_ercc),
        ('total UMIs',total_UMIs),
        ('detected genes',len(detected_genes))
    ])
    write_metrics(metricfile,metric_dict,metric_dict.keys())
    logger.info('Finished UMI counting and writing to disk')

def count_umis(gene_hash,primer_bed,tagged_bam,outfile_primer,outfile_gene,metricfile,logfile,cores):
    ''' Search for the design spe primers in the tagged bam
    and count molecular tags for each primer
    To do : Make this function shorter

    :param dict gene_hash: a dict of lists for storing gene annotations
    :param str primer_bed: a tsv file <chrom><start><stop><primer_seq><strand><gene>
    :param str tagged_bam: a UMI tagged bam file
    :param str outfile_primer: the output file for counts on primer level
    :param str outfile_gene: the output file for counts on a gene level
    :param str metricfile: file to write primer finding stats
    '''
    ## Set up logging
    logger = logging.getLogger("count_umis")
    logger.setLevel(logging.DEBUG)
    LOG = logging.FileHandler(logfile)
    logger.addHandler(LOG)
    ## Variable Initialization
    primer_info = defaultdict(list)
    primer_tree = defaultdict(lambda:IntervalTree())
    umi_counter = defaultdict(lambda:defaultdict(int))
    umi_counter_gene = defaultdict(lambda:defaultdict(int))
    patterns = defaultdict(lambda:defaultdict(list))
    ercc_used_unique=0
    ercc_used_multimapped=0
    multimapped = 0
    multimapped_used = 0
    primer_miss = 0
    primer_offtarget =0
    primer_mismatch = 0
    unmapped = 0
    endo_seq_miss=0
    endo_seq_miss_ercc=0
    num_reads_used_unique = 0    
    total_UMIs = 0
    ## Read the primer file and create an interval tree data structure
    with open(primer_bed) as IN:
        for line in IN:
            chrom,five_prime,three_prime,seq,strand,gene,ensembl_id = line.strip('\n').split('\t')
            sequence = Seq(seq)
            revcomp_seq = sequence.reverse_complement().tostring()
            if gene.startswith('ERCC-'): ## Update chrom
                chrom = gene
            if strand == '0':
                strand = '1'
                start = int(five_prime)
		stop = int(three_prime) + 1  ## Incrementing by 1 since interval tree assumes stop coordinate to be non-inclusive
                expression = regex.compile(r'^(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}[ACGTN]*$'%seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])
            else:
                strand = '-1'
                start = int(three_prime)
                stop = int(five_prime) + 1
                expression = regex.compile(r'^[ACGTN]*(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}$'%revcomp_seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])
            primer_info[seq] = [ensembl_id,seq,gene,strand,chrom,five_prime,three_prime]

    ## Iterate over the bam in chunks and process the results in parallel
    ## The chunking here is mainly to stay within memory bound for very large bam files
    p = Pool(cores)
    i = 1
    func = partial(find_primer,primer_tree)
    for chunks in iterate_bam_chunks(tagged_bam,chunks=10000000):
        find_primer_results = p.map(func,chunks)
        for info in find_primer_results:
            primer,umi,count,nh = info
            if nh>1:
                multimapped+=1
            if count == 0:
                if primer == 'Unknown_Chrom':
                    primer_offtarget+=1
                elif primer == 'Unmapped':
                    unmapped+=1
                elif primer == 'Unknown_Regex':
                    primer_mismatch+=1
                elif primer == 'Unknown_Loci':
                    primer_miss+=1
                else:
                    gene = primer_info[primer][2]
                    if gene.startswith('ERCC-'):
                        endo_seq_miss_ercc+=1
                    endo_seq_miss+=1
            else:
                gene = primer_info[primer][2]
                if nh > 1:
                    if gene.startswith('ERCC-'):                        
                        ercc_used_multimapped+=1
                    multimapped_used+=1
                else:
                    num_reads_used_unique+=1
                    if gene.startswith('ERCC-'):
                        ercc_used_unique+=1
                umi_counter[primer][umi]+=1
                umi_counter_gene[gene][umi]+=1
        i+=1
    p.close()
    p.join()
    ## Print output results
    seen = []
    detected_genes=0
    with open(outfile_primer,'w') as OUT1,open(outfile_gene,'w') as OUT2 :
        for primer in primer_info:
            ensembl_id,seq,gene,strand,chrom,five_prime,three_prime = primer_info[primer]
            if primer in umi_counter:
                umi_count = len(umi_counter[primer])
            else:
                umi_count = 0
            OUT1.write(ensembl_id+"\t"+gene+"\t"+strand+"\t"+str(chrom)+"\t"+str(five_prime)+"\t"+str(three_prime)+"\t"+seq+"\t"+str(umi_count)+"\n")
            if gene not in seen: ## Genes will be repeated for multiple primers, since results are already accumulated , only write once for a gene
                if gene in umi_counter_gene:
                    umi_count_gene = len(umi_counter_gene[gene])
                else:
                    umi_count_gene = 0
                total_UMIs+=umi_count_gene                
                OUT2.write(ensembl_id+"\t"+gene+"\t"+strand+"\t"+str(chrom)+"\t"+str(five_prime)+"\t"+str(three_prime)+"\t"+str(umi_count_gene)+"\n")
                if not gene.startswith('ERCC') and umi_count_gene > 0:
                    detected_genes+=1
                seen.append(gene)
                
    primers_found = len(umi_counter)    
    ## Write metrics
    num_reads_mapped_ercc = endo_seq_miss_ercc + ercc_used_unique + ercc_used_multimapped
    num_reads_mapped_genome = num_reads_used_unique + multimapped_used + endo_seq_miss + \
                              primer_mismatch + primer_miss + \
                              primer_offtarget + endo_seq_miss_ercc - \
                              num_reads_mapped_ercc
    
    num_reads_used_genome_unique = num_reads_used_unique - ercc_used_unique
    num_reads_used_genome_multimapped = multimapped_used - ercc_used_multimapped

    metric_dict = OrderedDict([
        ('reads dropped, not mapped to genome',unmapped),
        ('reads dropped, off target',primer_offtarget+primer_miss),
        ('reads dropped, primer not identified at read start',primer_mismatch),
        ('reads dropped, less than 25 bp endogenous seq after primer',endo_seq_miss),
        ('reads used, aligned to genome, multiple loci',num_reads_used_genome_multimapped),
        ('reads used, aligned to genome, unique loci',num_reads_used_genome_unique),
        ('reads used, aligned to ERCC, multiple loci',ercc_used_multimapped),
        ('reads used, aligned to ERCC, unique loci',ercc_used_unique),        
        ('detected genes',detected_genes),
        ('total UMIs',total_UMIs)
        ])
    
    write_metrics(metricfile,metric_dict,metric_dict.keys())
