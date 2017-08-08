from intervaltree import IntervalTree
import itertools
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
        yield [(read.seq, read.is_reverse, read.alen, chroms[read.tid]['SN'], read.pos, read.cigarstring, read.get_tag('mi')) for read in reads]

def count_mts(primer_bed,tagged_bam,outfile,metricfile):
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

    p = Pool(8)
    i = 1
    func = partial(find_primer,primer_tree)
    ## Iterate over the bam in chunks and process the results in parallel
    ## The chunking here is mainly to stay within memory bound for very large bam files
    ## The default chunk size is : 750,000 reads
    for chunks in iterate_bam_chunks(tagged_bam,chunks=1000000):
        print "Processing chunk : %i"%i
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

    ## Write metrics
    metric_dict = {
        'num_reads_primer_found':primer_found,
        'num_reads_primer_offtarget':primer_offtarget,
        'num_reads_primer_mismatch':primer_mismatch,
        'num_reads_primer_off_loci':primer_miss,
        'num_reads_unmapped':unmapped,
        'num_reads_endogenous_seq_not_matched':endo_seq_miss
    }
    write_metrics(metricfile,metric_dict)
    ## Print output results
    with open(outfile,'w') as OUT:
        for primer in primer_info:
            chrom,start,stop,seq,revcomp,strand,gene = primer_info[primer]
            mt_count = len(mt_counter[primer])
            print >> OUT,chrom+'\t'+start+'\t'+stop+'\t'+seq+'\t'+strand+'\t'+gene+'\t'+str(mt_count)

