import itertools
import sys
from guppy import hpy
import pysam
from pathos.multiprocessing import ProcessingPool as Pool
import regex
from collections import defaultdict,OrderedDict
from functools import Partial
from Bio.Seq import Seq

## Modules from this project
from find_primer import find_primer

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
    for reads in grouper(IN.fetch(),chunks):
        yield [(read.seq, chroms[read.tid]['SN'], read.get_tag('mi')) for read in reads]

def count_mts(primer_bed,tagged_bam,outfile):
    ''' Search for the design spe primers in the tagged bam
    and count molecular tags for each primer
    To do : Make this function shorter

    :param str primer_bed: a tsv file <chrom><start><stop><primer_seq><strand><gene>
    :param str tagged_bam: a UMI tagged bam file
    :param str outfile: the output file
    '''
    #hp = hpy()
    #print "Heap at the begining of the function\n", hp.heap()

    primer_dict = defaultdict(lambda:defaultdict(list))
    mt_counter = defaultdict(lambda:defaultdict(int))
    patterns = defaultdict(lambda:defaultdict(list))
    miss = 0

    with open(primer_bed) as IN:
        for line in IN:
            chrom,start,stop,seq,strand,gene = line.strip('\n').split('\t')
            sequence = Seq(seq)
            revcomp_seq = sequence.reverse_complement().tostring()
            primer_dict[chrom][seq] = [chrom,start,stop,seq,revcomp_seq,strand,gene]
            if strand == '+':
                #d<=1,i<=1,s<=3,3d+3i+1s<=7  A more complicated regex for future reference
                patterns[chrom][seq] = regex.compile(r'^[ACGTN]{0,5}(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}[ACGTN]*'%seq)
                #patterns[chrom][seq] = regex.compile(r'(%s)[ACGTN]*'%seq)
            else:
                #patterns[chrom][seq] = regex.compile(r'[ACGTN]*(%s)'%revcomp_seq)
                patterns[chrom][seq] = regex.compile(r'^[ACGTN]*(%s){d<=2,i<=2,s<=2,1d+1i+1s<=3}[ACGTN]{0,5}$'%revcomp_seq)

    p = Pool(8)
    i = 1
    #fuzzyness = tre.Fuzzyness(delcost=1,inscost=1,maxcost=3,subcost=1, maxdel=2,maxerr=3,maxins=2,maxsub=2)
    func = partial(find_primer,primer_dict,patterns)
    ## Iterate over the bam in chunks and process the results in parallel
    ## The chunking here is mainly to stay within memory bound for very large bam files
    ## The default chunk size is : 750,000 reads
    for chunks in iterate_bam_chunks(tagged_bam,chunks=1000000):
        print "Processing chunk : %i"%i
        find_primer_results = p.map(func,chunks)
        for info in find_primer_results:
            primer,mt,count = info
            if count == 0:
                miss+=1
            else:
                mt_counter[primer][mt]+=1
        i+=1
    p.close()
    p.join()
    print "Num reads not matched : %s"%miss

    ## Print output results
    with open(outfile,'w') as OUT:
        for chrom in primer_dict:
            for primer in primer_dict[chrom]:
                chrom,start,stop,seq,revcomp,strand,gene = primer_dict[chrom][primer]
                mt_count = len(mt_counter[primer])
                print >> OUT,chrom+'\t'+start+'\t'+stop+'\t'+seq+'\t'+strand+'\t'+gene+'\t'+str(mt_count)
    #print "Heap at the end of the function\n", hp.heap()
