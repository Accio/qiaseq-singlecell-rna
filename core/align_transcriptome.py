import subprocess
import sys
import pysam
import editdistance
import regex
from collections import defaultdict,OrderedDict

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

def find_primer(primer_dict,read_sequence,patterns):
    '''
    '''

    for primer in primer_dict:
        if regex.match(patterns[primer],read_sequence):
            return primer_dict[primer]
    return -1

def count_mts(primer_bed,tagged_bam,outfile):
    ''' Search for the primers in the bedfile in the tagged bam
    and count mt for each primer
    :param str primer_bed: a tsv file <chrom><start><stop><primer_seq><strand><gene>
    :param str tagged_bam: a UMI tagged bam file
    :param str outfile: the output file
    '''

    primer_dict = OrderedDict()
    mt_counter = defaultdict(lambda:defaultdict(int))
    miss = 0
    pattern = ('^(%s){e<=2}[ACGTN]*')
    patterns = {}
    with open(primer_bed) as IN:
        for line in IN:
            chrom,start,stop,seq,strand,gene = line.strip('\n').split('\t')
            primer_dict[seq] = [chrom,start,stop,seq,gene]
            patterns[seq] = regex.compile(pattern%seq)

    with pysam.AlignmentFile(tagged_bam,'rb') as IN:
        for read in IN:
            mt = read.get_tag('mi')
            match = find_primer(primer_dict,read.seq,patterns)
            if match == -1:
                miss+=1
            else:
                primer_seq = match[-2]
                mt_counter[primer_seq][mt]+=1

    print "Num reads not matched : %s"%miss

    with open(outfile,'w') as OUT:
        for key in primer_dict:
            chrom,start,stop,seq,gene = primer_dict[key]
            mt_count = len(mt_counter[key].keys())
            print >> OUT,chrom+'\t'+start+'\t'+stop+'\t'+seq+'\t'+gene+'\t'+str(mt_count)


if __name__ == '__main__':
    #star_alignment("/qgen/home/jdeng/download/STAR-2.5.0b/bin/Linux_x86_64_static/STAR",sys.argv[1],sys.argv[2],"--runMode alignReads --genomeLoad LoadAndKeep --runThreadN 6 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --alignIntronMax 200000 --alignMatesGapMax 200000 --alignSJDBoverhangMin 16 --sjdbOverhang 149 --outSAMunmapped Within --outSAMprimaryFlag AllBestScore --outSAMmultNmax 1",sys.argv[3])
    #annotate_bam_umis(sys.argv[1],sys.argv[2],sys.argv[3])
    count_mts(sys.argv[1],sys.argv[2],sys.argv[3])
