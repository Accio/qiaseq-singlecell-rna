from intervaltree import IntervalTree
import itertools
import logging
import sys
from guppy import hpy
import pysam
from pathos.multiprocessing import ProcessingPool as Pool
import regex
from collections import defaultdict,OrderedDict
from functools import partial
from Bio.Seq import Seq

## Modules from this project
from find_primer import find_primer
from demultiplex_cells import write_metrics
from extract_multiplex_region import open_by_magic

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

    :param str tagged_bam: the input bam file with MT tags
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

def create_gene_tree(annotation_gtf,merge_coordinates=False):
    '''
    :param str annotation_gtf : a gtf file for identifying genic regions
    '''
    gene_tree = defaultdict(lambda:defaultdict(IntervalTree))
    genes = defaultdict(list)
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                start = int(contents[3])
                end = int(contents[4])
                strand = contents[6]
                info = contents[-1]
                gene = info.split(';')[3].split()[1].strip('\"')
                gene_type = info.split(';')[1].split()[1].strip('\"')

                if gene == None or gene_type == None:
                    raise Exception(
                        "Failed Parsing annotation file :{annotation}".format(
                            annotation=annotation_gtf))

                if merge_coordinates: ## Create a coordinate set which is merged to include the largest interval possible
                    if gene in genes: ## Gene has been seen before
                        if genes[gene][2] != chrom: ## Different chromosome
                            genes[gene] = (start,end,chrom,strand,gene,gene_type)
                        else:
                            if start <= genes[gene][0]: ## Gene has unique start bases to add
                                if end >= genes[gene][1]: ## Gene has unique end bases to add
                                    genes[gene] = (start,end,chrom,strand,gene,gene_type)
                                else:
                                    if end < genes[gene][0]: ## Check if this interval is dijoint from the previous one
                                        genes[gene].append((start,end,chrom,strand,gene,gene_type)) ## Add a new interval
                                    else:
                                        genes[gene] = (start,genes[gene][1],chrom,strand,gene,gene_type) ## Merge intervals

                            else: ## Start is already spanned , check end
                                if end >= genes[gene][1]: ## Update end base position
                                    if start > genes[gene][1]: ## Start is greater than previously encountered gene's end
                                        genes[gene].append((start,end,chrom,strand,gene,gene_type)) ## Add a new interval
                                    else:
                                        genes[gene] = (genes[gene][0],end,chrom,strand,gene,gene_type) ## Merge intervals
                                else: ## No need to update anything
                                    continue
                    else:
                        genes[gene] = (start,end,chrom,strand,gene,gene_type)

                else: ## Store all intervals without mergeing
                    genes[gene].append((start,end,chrom,strand,gene,gene_type))

    ## Build a primer tree to store gene info
    for gene in genes:
        for info in genes[gene]:
            start,end,chrom,strand,gene,gene_type = info
            gene_tree[chrom][strand].addi(start,end+1,info)

    print "Interval tree created with {ngenes}".format(ngenes=len(genes))
    return gene_tree

def overlap(x1,x2,y1,y2):
    ''' Compute overlap between two coordinates
    '''
    return max(0,min(x2,y2) - max(x1,y1))

def find_gene(gene_tree,read_tup):
    ''' Annotate the given read with a gene

    :param interval treee gene_tree: efficient interval tree data structure storing coordinates and gene information
    :param tuple read_tup: a tuple of read information
    :return a tuple containing gene and mt info
    :rtype tuple
    '''

    overlap_threshold = 0
    logger = logging.getLogger("count_umis")
    read_id,read_sequence,read_is_reverse,read_len,read_chrom,read_pos,read_cigar,mt,nh = read_tup
    if read_chrom == "*":
        logger.info("{read_id}: Unmapped".format(read_id=read_id))
        return ('Unmapped',mt,0)
    if nh != 1: ## Multimapped read
        logger.info("{read_id}: Multimapped".format(read_id=read_id))
        return ('Multimapped',mt,0)
    if 'ERCC' in read_chrom:
        logger.info("{read_id}: Mapped to {ercc}".format(read_id=read_id,ercc=read_chrom))
        return (read_chrom,mt,0)
    if read_chrom not in gene_tree:
        logger.info("{read_id}: Chromosome {chrom} was not present in annotation gene interval".format(read_id=read_id,chrom=read_chrom))
        return ('Unknown_Chrom',mt,0)

    ## Search the interval tree
    read_end = return_read_end_pos(read_pos,read_cigar)
    if read_is_reverse:
        res = gene_tree[read_chrom]['-'].search(read_pos,read_end)
    else:
        res = gene_tree[read_chrom]['+'].search(read_pos,read_end)

    if res: ## If the search was successful
        read_end = return_read_end_pos(read_pos,read_cigar,flag=True)
        num_hits = len(res)
        logger.info("{}".format(num_hits))
        if num_hits > 1:
            logger.info("{read_id}: intersected with {hits} genic regions".format(read_id=read_id,hits=num_hits))
            ## Choose the closest 3' location gene
            prev = None
            for result in res:
                if not prev:
                    ## Check overlap
                    prev_o = float(overlap(read_pos,read_end,result.data[0],result.data[1]))/read_len
                    if prev_o > overlap_threshold:
                        prev = result.data
                        logger.info("{read_id}: Picked {gene} as default".format(read_id=read_id,gene=prev[4]))
                    else:
                        logger.info("{read_id}: {gene} had {overlap} overlap, failed overlap criteria".format(read_id=read_id,overlap=prev_o,gene=result.data[4]))
                else:
                    o = float(overlap(read_pos,read_end,result.data[0],result.data[1]))/read_len
                    if o < overlap_threshold:
                        continue
                    if read_is_reverse:
                        three_prime_prev = prev[0]
                        three_prime_new = result.data[0]
                        three_prime_read = read_pos
                    else:
                        three_prime_prev = prev[1]
                        three_prime_new = result.data[1]
                        three_prime_read = read_end

                    gene_type_new = result.data[-1]
                    gene_type_old = prev[-1]
                    diff_three_prime_prev = abs(three_prime_prev - three_prime_read)
                    diff_three_prime_new = abs(three_prime_new - three_prime_read)

                    if gene_type_new != gene_type_old:
                        if gene_type_new != 'protein_coding':
                            if gene_type_old == 'protein_coding':
                                logger.info("{read_id}: Checking {gene} , not chosen as it is not protein_coding".format(read_id=read_id,gene=result.data[4]))
                                continue
                        else:
                            logger.info("{read_id}: Chose {gene2} over {gene1} because it is protein coding".format(read_id=read_id,gene1=prev[4],gene2=result.data[4]))
                            prev = result.data
                            prev_o = o
                            continue

                    if diff_three_prime_prev == diff_three_prime_new:
                        if o < prev_o: ## Look at overlaps
                            logger.info("{read_id}: Picked {gene1}: 3_prime_diff={diff1} over {gene2}: 3_prime_diff={diff2} because of greater overlap".format(read_id=read_id,gene1=result.data[4],diff1=diff_three_prime_new,gene2=prev[4],diff2=diff_three_prime_prev))
                            prev_o = o
                            prev = result.data
                    elif diff_three_prime_prev > diff_three_prime_new:
                        logger.info("{read_id}: Picked {gene1}: 3_prime_diff={diff1} over {gene2}: 3_prime_diff={diff2}".format(read_id=read_id,gene1=result.data[4],diff1=diff_three_prime_new,gene2=prev[4],diff2=diff_three_prime_prev))
                        prev = result.data
                        prev_o = o
            if prev:
                logger.info("{read_id}: Picked {gene1}".format(read_id=read_id,gene1=prev[4]))
                return (prev,mt,1)
            else:
                logger.info("{read_id}: No genes matched overlap criteria".format(read_id=read_id))
                return ('Unknown',mt,0)
        else:
            result = res.pop()
            logger.info("{read_id}: intersected with {gene} only".format(read_id=read_id,gene=result.data[4]))
            o = float(overlap(read_pos,read_end,result.data[0],result.data[1]))/read_len
            if o < overlap_threshold:
                logger.info("{read_id}: {gene} failed overlap criteria".format(read_id=read_id,gene=result.data[4]))
                return ('Unknown',mt,0)
            else:
                return (result.data,mt,1)
    else: ## Could not find loci in gene tree
        logger.info("{read_id} was not found in the annotation gene interval tree".format(read_id=read_id))
        return ('Unknown',mt,0)

def return_read_end_pos(read_pos,cigar,flag=False,pattern=regex.compile('([0-9]+)([A-Z])'),cigars_to_ignore=['I','S','H','P']):
    ''' Return the end of the read
    '''
    expanded_cigar = []
    bases=0
    if flag:
        bases_to_ignore = ['S','I','H','P','N']
    for num_bases,cigar_char in regex.findall(pattern,cigar):
        if cigar_char in cigars_to_ignore:
            continue
    else:
        bases+=int(num_bases)
    return read_pos+bases

def count_umis_wts(annotation_gtf,tagged_bam,outfile,metricfile,logfile,cores=3):
    ''' Count MTs for each gene in the input the tagged_bam file

    :param str annotation_gtf : a gtf file for identifying genic regions
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

    found = 0
    miss_chr = 0
    unmapped = 0
    not_annotated = 0
    multimapped = 0
    ercc=0
    ercc_names = ['ERCC-00002','ERCC-00003','ERCC-00004','ERCC-00009','ERCC-00012','ERCC-00013','ERCC-00014','ERCC-00016','ERCC-00017','ERCC-00019','ERCC-00022','ERCC-00024','ERCC-00025','ERCC-00028','ERCC-00031','ERCC-00033','ERCC-00034','ERCC-00035','ERCC-00039','ERCC-00040','ERCC-00041','ERCC-00042','ERCC-00043','ERCC-00044','ERCC-00046','ERCC-00048','ERCC-00051','ERCC-00053','ERCC-00054','ERCC-00057','ERCC-00058','ERCC-00059','ERCC-00060','ERCC-00061','ERCC-00062','ERCC-00067','ERCC-00069','ERCC-00071','ERCC-00073','ERCC-00074','ERCC-00075','ERCC-00076','ERCC-00077','ERCC-00078','ERCC-00079','ERCC-00081','ERCC-00083','ERCC-00084','ERCC-00085','ERCC-00086','ERCC-00092','ERCC-00095','ERCC-00096','ERCC-00097','ERCC-00098','ERCC-00099','ERCC-00104','ERCC-00108','ERCC-00109','ERCC-00111','ERCC-00112','ERCC-00113','ERCC-00116','ERCC-00117','ERCC-00120','ERCC-00123','ERCC-00126','ERCC-00130','ERCC-00131','ERCC-00134','ERCC-00136','ERCC-00137','ERCC-00138','ERCC-00142','ERCC-00143','ERCC-00144','ERCC-00145','ERCC-00147','ERCC-00148','ERCC-00150','ERCC-00154','ERCC-00156','ERCC-00157','ERCC-00158','ERCC-00160','ERCC-00162','ERCC-00163','ERCC-00164','ERCC-00165','ERCC-00168','ERCC-00170','ERCC-00171']
    ercc_info = {}
    N=10000000

    ## Build an efficient data structure to query gene annotations
    gene_tree = create_gene_tree(annotation_gtf)
    ## Store mt counts for each gene
    mt_counter = defaultdict(lambda:defaultdict(int))
    logger.info('Using {} cores'.format(cores))
    p = Pool(cores)
    func = partial(find_gene,gene_tree)
    for reads in iterate_bam_chunks(tagged_bam,chunks=N):
        logger.info('Reading {} reads in memory to find genes'.format(N))
        find_gene_results = p.map(func,reads)
        for info in find_gene_results:
            gene_info,mt,count = info
            if count == 0:
                if gene_info == 'Unknown_Chrom':
                    miss_chr+=1
                elif gene_info == 'Unmapped':
                    unmapped+=1
                elif gene_info == 'Multimapped':
                    multimapped+=1
                elif gene_info == 'Unknown':
                    not_annotated+=1
                else:
                    ercc+=1
                    info = ('N/A','N/A','N/A','N/A',gene_info,'N/A')
                    ercc_info[gene_info] = info
                    mt_counter[info][mt]+=1
            else:
                mt_counter[gene_info][mt]+=1
                found+=1
    p.close()
    p.join()
    logger.info('Finished gene finding. Writing metrics and counting MTs')
    ## Write metrics
    metric_dict = OrderedDict([
        ('num_reads_mapped',found+miss_chr+not_annotated+multimapped),
        ('num_reads_not_annotated',not_annotated),
        ('num_reads_unknown_chrom', miss_chr),
        ('num_reads_multimapped',multimapped),
        ('num_reads_used',found),
        ('num_reads_used_ercc',ercc),
        ('num_genes_annotated',len(mt_counter))
    ])
    write_metrics(metricfile,metric_dict)
    ## Print output results
    ## Write gene counts
    with open(outfile,'w') as OUT:
        for chrom in gene_tree:
            for strand in gene_tree[chrom]:
                for gene_info in gene_tree[chrom][strand]:
                    start,end,chrom,strand,gene,gene_type = gene_info.data
                    #loci = chrom+':'+start+'-'+end
                    if gene_info.data in mt_counter:
                        mt_count = len(mt_counter[gene_info.data])
                    else:
                        mt_count = 0
                    print >> OUT,chrom+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+gene+\
                        '\t'+gene_type+'\t'+str(mt_count)
        ## Write ERCC counts
        for name in ercc_names:
            if name in ercc_info:                
                gene_info = ercc_info[name]
                start,end,chrom,strand,gene,gene_type = gene_info
                mt_count = len(mt_counter[gene_info])
            else:
                start=end=chrom=strand=gene_type='N/A'
                gene = name
                mt_count = 0
            print >> OUT,chrom+'\t'+str(start)+'\t'+str(end)+'\t'+strand+'\t'+gene+\
                '\t'+gene_type+'\t'+str(mt_count)
    logger.info('Finished MT counting and writing to disk')

def count_umis(primer_bed,tagged_bam,outfile,metricfile):
    ''' Search for the design spe primers in the tagged bam
    and count molecular tags for each primer
    To do : Make this function shorter

    :param str primer_bed: a tsv file <chrom><start><stop><primer_seq><strand><gene>
    :param str tagged_bam: a UMI tagged bam file
    :param str outfile: the output file
    :param str metricfile: file to write primer finding stats
    '''
    primer_info = defaultdict(list)
    primer_tree = defaultdict(lambda:IntervalTree())
    #primer_dict = defaultdict(lambda:defaultdict(list))
    mt_counter = defaultdict(lambda:defaultdict(int))
    patterns = defaultdict(lambda:defaultdict(list))
    primer_miss=0
    primer_offtarget=0
    primer_mismatch=0
    unmapped=0
    endo_seq_miss=0
    primer_found = 0
    with open(primer_bed) as IN:
        for line in IN:
            chrom,start,stop,seq,strand,gene = line.strip('\n').split('\t')
            sequence = Seq(seq)
            revcomp_seq = sequence.reverse_complement().tostring()
            #primer_dict[chrom][seq] = [chrom,start,stop,seq,revcomp_seq,strand,gene]
            primer_info[seq] = [chrom,start,stop,seq,revcomp_seq,strand,gene]
            if strand == '+':
                expression = regex.compile(r'^(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}[ACGTN]*$'%seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])
            else:
                expression = regex.compile(r'^[ACGTN]*(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}$'%revcomp_seq)
                primer_tree[chrom].addi(int(start),int(stop),[expression,seq])

    p = Pool(4)
    i = 1
    func = partial(find_primer,primer_tree)
    ## Iterate over the bam in chunks and process the results in parallel
    ## The chunking here is mainly to stay within memory bound for very large bam files
    ## The default chunk size is : 750,000 reads
    for chunks in iterate_bam_chunks(tagged_bam,chunks=1000000):
        #print "Processing chunk : %i"%i
        find_primer_results = p.map(func,chunks)
        for info in find_primer_results:
            primer,mt,count = info
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
                    endo_seq_miss+=1
            else:
                primer_found+=1
                mt_counter[primer][mt]+=1
        i+=1
    p.close()
    p.join()

    primers_found = len(mt_counter)
    ## Write metrics
    metric_dict = OrderedDict([
        ('num_primers_found',primers_found),
        ('num_reads_primer_offtarget',primer_offtarget),
        ('num_reads_primer_mismatch',primer_mismatch,)
        ('num_reads_primer_off_loci',primer_miss,)
        ('num_reads_unmapped',unmapped),
        ('num_reads_endogenous_seq_not_matched',endo_seq_miss),
        ('num_reads_primer_found',primer_found),
    ])
    write_metrics(metricfile,metric_dict)

    ## Print output results
    with open(outfile,'w') as OUT:
        for primer in primer_info:
            chrom,start,stop,seq,revcomp,strand,gene = primer_info[primer]
            mt_count = len(mt_counter[primer])
            print >> OUT,chrom+'\t'+start+'\t'+stop+'\t'+strand+'\t'+gene+'\t'+seq+'\t'+str(mt_count)

