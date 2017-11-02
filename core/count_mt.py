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
from extract_multiplex_region import open_by_magic
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
            to_yield.append((read.qname,read.seq, read.is_reverse, read.alen,chromosome,
                             read.pos, read.cigarstring, read.get_tag('mi'),
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
    miss_chr = 0
    unmapped = 0
    not_annotated = 0
    multimapped = 0
    ercc  = 0
    ercc_multimapped = 0
    total_UMIs = 0
    ## To do : Store an ERCC file to output correct coordinates
    ercc_names = ['ERCC-00002','ERCC-00003','ERCC-00004','ERCC-00009','ERCC-00012','ERCC-00013','ERCC-00014','ERCC-00016','ERCC-00017','ERCC-00019','ERCC-00022','ERCC-00024','ERCC-00025','ERCC-00028','ERCC-00031','ERCC-00033','ERCC-00034','ERCC-00035','ERCC-00039','ERCC-00040','ERCC-00041','ERCC-00042','ERCC-00043','ERCC-00044','ERCC-00046','ERCC-00048','ERCC-00051','ERCC-00053','ERCC-00054','ERCC-00057','ERCC-00058','ERCC-00059','ERCC-00060','ERCC-00061','ERCC-00062','ERCC-00067','ERCC-00069','ERCC-00071','ERCC-00073','ERCC-00074','ERCC-00075','ERCC-00076','ERCC-00077','ERCC-00078','ERCC-00079','ERCC-00081','ERCC-00083','ERCC-00084','ERCC-00085','ERCC-00086','ERCC-00092','ERCC-00095','ERCC-00096','ERCC-00097','ERCC-00098','ERCC-00099','ERCC-00104','ERCC-00108','ERCC-00109','ERCC-00111','ERCC-00112','ERCC-00113','ERCC-00116','ERCC-00117','ERCC-00120','ERCC-00123','ERCC-00126','ERCC-00130','ERCC-00131','ERCC-00134','ERCC-00136','ERCC-00137','ERCC-00138','ERCC-00142','ERCC-00143','ERCC-00144','ERCC-00145','ERCC-00147','ERCC-00148','ERCC-00150','ERCC-00154','ERCC-00156','ERCC-00157','ERCC-00158','ERCC-00160','ERCC-00162','ERCC-00163','ERCC-00164','ERCC-00165','ERCC-00168','ERCC-00170','ERCC-00171']
    ercc_info = {}
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
                else:
                    ercc+=1
                    if nh>1:
                        ercc_multimapped+=1
                        multimapped+=1
                    else:
                        info = ('N/A','N/A','N/A','N/A',gene_info,'N/A')
                        ercc_info[gene_info] = info
                        umi_counter[info][umi]+=1
                    found+=1
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
    with open(outfile,'w') as OUT:
        for chrom in gene_tree:
            for strand in gene_tree[chrom]:
                for gene_info in gene_tree[chrom][strand]:
                    start,end,chrom,strand,gene,gene_type = gene_info.data
                    if gene_info.data in umi_counter:
                        umi_count = len(umi_counter[gene_info.data])
                    else:
                        umi_count = 0
                    total_UMIs+=umi_count
                    print >> OUT,chrom+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+gene+\
                        '\t'+gene_type+'\t'+str(umi_count)
        ## Write ERCC counts
        for name in ercc_names:
            if name in ercc_info:                
                gene_info = ercc_info[name]
                start,end,chrom,strand,gene,gene_type = gene_info
                umi_count = len(umi_counter[gene_info])
            else:
                start=end=chrom=strand=gene_type='N/A'
                gene = name
                umi_count = 0
            total_UMIs+=umi_count
            print >> OUT,chrom+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+gene+\
                '\t'+gene_type+'\t'+str(umi_count)
    ## Write metrics
    metric_dict = OrderedDict([
        ('num_reads_mapped',found+miss_chr+not_annotated),
        ('num_reads_mapped_ercc',ercc),
        ('num_reads_not_annotated',not_annotated),
        ('num_reads_unknown_chrom', miss_chr),
        ('num_reads_multimapped',multimapped),
        ('num_reads_uniquely_mapped',found+miss_chr+not_annotated-multimapped),
        ('num_reads_used',found-multimapped),
        ('num_reads_used_ercc',ercc-ercc_multimapped),
        ('num_umis_used',total_UMIs),
        ('num_genes_annotated',len(umi_counter))
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
    ercc = 0
    ercc_used=0
    multimapped = 0
    primer_miss = 0
    primer_offtarget =0
    primer_mismatch = 0
    unmapped = 0
    endo_seq_miss=0
    num_reads_used = 0    
    total_UMIs = 0
    ## Read the primer file and create an interval tree data structure
    with open(primer_bed) as IN:
        for line in IN:
            chrom,five_prime,three_prime,seq,strand,gene = line.strip('\n').split('\t')
            sequence = Seq(seq)
            revcomp_seq = sequence.reverse_complement().tostring()
            if gene.startswith('ERCC-'): ## Update chrom
                chrom = gene
            if strand == '0':
                strand = '+'
                start = int(five_prime)
		stop = int(three_prime) + 1  ## Incrementing by 1 since interval tree assumes stop coordinate to be non-inclusive
                expression = regex.compile(r'^(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}[ACGTN]*$'%seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])
            else:
                strand = '-'
                start = int(three_prime)
                stop = int(five_prime) + 1
                expression = regex.compile(r'^[ACGTN]*(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}$'%revcomp_seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])
            primer_info[seq] = [chrom,start,stop,seq,revcomp_seq,strand,gene]

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
                    gene = primer_info[primer][-1]
                    if gene.startswith('ERCC-'):
                        ercc+=1
                    endo_seq_miss+=1
            else:
                num_reads_used+=1
                umi_counter[primer][umi]+=1
                gene = primer_info[primer][-1]
                umi_counter_gene[gene][umi]+=1
                if gene.startswith('ERCC-'):
                    ercc_used+=1
        i+=1
    p.close()
    p.join()
    ## Print output results
    with open(outfile_primer,'w') as OUT1,open(outfile_gene,'w') as OUT2 :
        for primer in primer_info:
            chrom,start,stop,seq,revcomp,strand,gene = primer_info[primer]
            if primer in umi_counter:
                umi_count = len(umi_counter[primer])
            else:
                umi_count = 0
            OUT1.write(chrom+'\t'+str(start)+'\t'+str(stop)+'\t'+strand+'\t'+gene+'\t'+seq+'\t'+str(umi_count)+'\n')
            if gene in umi_counter_gene:
                umi_count = len(umi_counter_gene[gene])
            else:
                umi_count = 0
            total_UMIs+=umi_count    
            if len(gene_hash[gene]) != 6: ## Temp fix for dealing with edge cases where gene is not present in annotation file
                temp = ['N/A','N/A','N/A','N/A',gene,'N/A']
                gene_info = '\t'.join(temp)
            else:
                gene_info = '\t'.join(gene_hash[gene])                
            OUT2.write(gene_info+'\t'+str(umi_count)+'\n')
                
    primers_found = len(umi_counter)
    genes_found = len(umi_counter_gene)
    ## Write metrics
    num_reads_mapped = num_reads_used + endo_seq_miss + \
                       primer_mismatch + primer_miss + \
                       primer_offtarget + ercc
    num_reads_mapped_ercc = ercc + ercc_used
    num_reads_used_ercc = ercc_used
    
    metric_dict = OrderedDict([
        ('num_primers_found',primers_found),
        ('num_genes_found',genes_found),
        ('num_reads_mapped',num_reads_mapped),
        ('num_reads_mapped_ercc',num_reads_mapped_ercc),
        ('num_reads_multimapped',multimapped),
        ('num_reads_primer_offtarget',primer_offtarget),
        ('num_reads_primer_mismatch',primer_mismatch),
        ('num_reads_primer_off_loci',primer_miss),
        ('num_reads_unmapped',unmapped),
        ('num_reads_endogenous_seq_not_matched',endo_seq_miss),
        ('num_reads_used',num_reads_used),
        ('num_umis_used',total_UMIs)
    ])
    write_metrics(metricfile,metric_dict,metric_dict.keys())

