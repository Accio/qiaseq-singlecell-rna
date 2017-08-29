import subprocess
import pysam
import sys

def run_cmd(cmd):
    ''' Run a shell command
    :param str cmd: the command to run
    '''
    p = subprocess.Popen(cmd,shell=True)
    p.wait()
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode,cmd)
    
def star_load_index(star,genome_dir,program_options):
    ''' Load star index
    :param str star: path to the star executable
    :param str genome_dir: path to the genomedir
    :param str program_options: options to use with star
    '''
    cmd = star + ' --genomeDir %s'%genome_dir + ' ' + program_options
    run_cmd(cmd)
    
def star_remove_index(star,genome_dir,program_options):
    ''' Remove star index
    '''
    cmd = star + ' --genomeDir %s'%genome_dir + ' ' + program_options
    run_cmd(cmd)
    
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
    run_cmd(cmd)
    
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

    try:
        with pysam.AlignmentFile(in_bam,'rb') as IN, pysam.AlignmentFile(out_bam,'wb',template=IN) as OUT:
            for read in IN:
                temp_tags = read.tags
                tag = tag_hash[read.qname]
                temp_tags.append((tag_name,tag))
                read.tags = tuple(temp_tags)
                OUT.write(read)
    except Exception as e:
        if e.message == "file header is empty (mode='rb') - is it SAM/BAM format?": ## An empty bam file
            print "Empty bam file : %s"%in_bam
            return
        else:
            raise Exception(e)
    pysam.index(out_bam)

def annotate_bam_umi_sort_merge(multiplex_file,in_bam,out_bam,tag_name="mi"):
    ''' Explore option with STAR to keep umi in read id and parse it out when iterating bam again
    '''
    sort_cmd = "sort -k{column,column} -"
    return 0
