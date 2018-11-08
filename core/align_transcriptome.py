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
    
def star_alignment(star,genome_dir,output_dir,logfile,program_options,r1,r2=None):
    '''
    Wrapper function to call STAR aligner with appropriate options

    :param str star: path to the star executable
    :param str genome_dir: path to the dir with genome and index files
    :param str output_dir: path to the output directory
    :param str logfile : log to redirect stderr and stdout to
    :param str program_options: options to use with star
    :param str r1: path to r1 fastq
    :param str r2: path to r2 fastq , default=None
    :return: does not return anything
    '''
    cmd = star + ' --genomeDir %s'%genome_dir + ' ' + program_options + ' --outFileNamePrefix %s'%output_dir + \
    ' --readFilesIn %s'%r1
    if r2: ## if using R2 reads as well
        cmd = cmd + ' ' + r2
    # redirect stderr and stdout to log
    cmd = cmd + ' ' + '> {log} 2>&1'.format(log=logfile)
    run_cmd(cmd)
