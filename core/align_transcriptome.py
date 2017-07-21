import subprocess
import itertools
import sys
from guppy import hpy
import pysam
from pathos.multiprocessing import ProcessingPool as Pool
import editdistance
import regex
from collections import defaultdict,OrderedDict
from functools import partial
from Bio.Seq import Seq

def star_alignment(star,genome_dir,output_dir,program_options,r1,r2=None):
    '''
    Wrapper function to call STAR aligner with appropriate options

    :param str star: path to the star executable
    :param str genome_dir: path to the dir with genome and index files
    :param str output_dir: path to the output directory
    :param str program_options: options to use with star
    :param str r1: path to r1 fastq
    :param str r2: path to r2 fastq , default=None
    :return: does not return anything
    '''

    cmd = star + ' --genomeDir %s'%genome_dir + ' ' + program_options + ' --outFileNamePrefix %s'%output_dir + \
    ' --readFilesIn %s'%r1
    if r2: ## if using R2 reads as well
        cmd = cmd + ' ' + r2

    p = subprocess.Popen(cmd,shell=True)
    p.wait()
    if p.returncode: ## Non zero returncode
        raise subprocess.CalledProcessError(p.returncode,cmd)


def annotate_bam_umi(multiplex_file,in_bam,out_bam,tag_name="mi"):
    ''' Annotate a bam file with UMIs, a new UMI tag will be created

    :param str multiplex_file: tsv file <read_id> <cell_index> <mt>
    :param str in_bam: an aligned bam file
    :param str out_bam: output bam file to write
    :param str tag_name: the name of the tag to use for the UMI
    :returns: nothing
    '''

    tag_hash = {}
    with open(multiplex_file,'r') as IN:
        for line in IN:
            contents = line.rstrip('\n').split('\t')
            tag_hash[contents[0].split()[0][1:]] = contents[2]

    with pysam.AlignmentFile(in_bam,'rb') as IN, pysam.AlignmentFile(out_bam,'wb',template=IN) as OUT:
        for read in IN:
            temp_tags = read.tags
            tag = tag_hash[read.qname]
            temp_tags.append((tag_name,tag))
            read.tags = tuple(temp_tags)
            OUT.write(read)


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

def find_primer(primer_dict,patterns,read_tup):
    '''
    Find whether a read matches one of the SPE primers used for the
    sequencing experiment

    :param dict of dict of lists primer_dict: Contains primer info
    :param a dict of dict  of lists patterns: Contains the compiled regular expression
                                              for the primer
    :param tuple read_tup: (read_sequence,chromosome,MT)
    :returns: a tuple containing the primer and mt info and whether it was a match
    :rtype: tuple
    '''

    read_sequence,read_chrom,mt = read_tup
    if read_chrom not in primer_dict:
        print "Invalid Chromosome for primer : %s"%(read_chrom)
        return ('Unknown','Unknown',0)

    for primer in primer_dict[read_chrom]:
        if regex.match(patterns[read_chrom][primer],read_sequence):
            return (primer , mt, 1)

    ## Fall back search on other chromosomes
    for chrom in primer_dict:
        if chrom == read_chrom:
            continue
        for primer in primer_dict[chrom]:
            if regex.match(patterns[chrom][primer],read_sequence):
                return (primer, mt, 1)

    return (primer,mt,0)

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
                patterns[chrom][seq] = regex.compile(r'^(%s){e<3}[ACGTN]*'%seq)
            else:
                patterns[chrom][seq] = regex.compile(r'^[ACGTN]*(%s){e<3}$'%revcomp_seq)

    p = Pool(12)
    i = 1
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

if __name__ == '__main__':
    ## Need to add command line arguments here
    #star_alignment("/qgen/home/jdeng/download/STAR-2.5.0b/bin/Linux_x86_64_static/STAR",sys.argv[1],sys.argv[2],"--runMode alignReads --genomeLoad LoadAndKeep --runThreadN 6 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --alignIntronMax 200000 --alignMatesGapMax 200000 --alignSJDBoverhangMin 16 --sjdbOverhang 149 --outSAMunmapped Within --outSAMprimaryFlag AllBestScore --outSAMmultNmax 1",sys.argv[3])
    #annotate_bam_umis(sys.argv[1],sys.argv[2],sys.argv[3])
    count_mts(sys.argv[1],sys.argv[2],sys.argv[3])
