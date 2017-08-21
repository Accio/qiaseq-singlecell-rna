import re
import argparse
import gzip
import multiprocessing
import time
import regex
import sys
import io
from functools import partial
from collections import Counter

"""
Given a Read sequence extract the cellindex and molecular tag region.
The read sequence in question here originates from the 2nd mate pair
of a single cell rna sequencing dataset.

R2 sequence is formatted as below :
<vector_sequence><cell_id><molecular_tag><isolator><TTTTTTTT..>

vector_sequence :  Constant nucleotide sequence for a given readset
cell_id : Fixed length nucleodtide sequence used to demultiplex reads to individual cells
molecular_tag : Fixed length nucleodtide sequence used to infer if a read originated from the same input molecule
isolator : Constant nucleotide sequence adjacent to the MT tag to identify the start of a polyT region

The goal is to extract the <cell_id><molecular_tag> region ; Any matching performed needs to account for some amount
of error in the vector, cell_id and mt regions
"""


def open_by_magic(filename):
    '''
    Adapted from : http://stackoverflow.com/questions/18367511/how-do-i-automatically-handle-decompression-when-reading-a-file-in-python
    with modifications
    Uses the initial bytes of a file to detect the file compression.

    :param str filename: path to the input file
    :return: the appropriate file handle for reading
    :rtype: file object
    '''

    ## Add more magic strs here for various compressions
    magic_dict = {"\x1f\x8b\x08":gzip.open}
    max_len = max(len(x) for x in magic_dict)
    with open(filename) as f:
        file_start = f.read(max_len)
        for magic,fn in magic_dict.items():
            if file_start.startswith(magic):
                return io.BufferedReader(fn(filename))
            return open(filename,'r') ## Otherwise just a regular file

def iterate_fastq(read2_fastq):
    '''
    Iterate and parse a fastq file

    :param str read2_fastq: the 2nd mate pair fastq file location
    :return: a (read_id,sequence) generator
    :rtype: generator of tuples
    '''

    with open_by_magic(read2_fastq) as IN:
        while True:
            yield (IN.next().rstrip('\n'),IN.next().rstrip('\n'))
            IN.next()
            IN.next()

def find_motif(motif,read_tup):
    '''
    Function to extract the <cell_index><mt> region

    :param str motif: regular expression to match read sequence
    :param tuple read_tup : (read_id,read_seq) tuple
    :return: multiplex region matching the motif
    :rtype: str
    '''

    read_id,read_seq = read_tup
    match = regex.match(motif,read_seq)
    if match:
        return (read_id,match.group(3))
    else:
        return (read_id,None)

def print_result(regions,outfile,cell_index_len):
    '''
    Output the results as a tsv file

    :param list of tuples regions: tuple of (read_id,'<cell_index><mt>')
    :param str outfile: path to the outputfile
    :return: nothing
    :rtype:
    '''
    i = 0
    j = 0
    with open(outfile,'w') as OUT:
        OUT.write("read_id\tcell_index\tmt\n")
        for read_id,multiplex_region in regions:
            if multiplex_region:
                i+=1
                cell_index = multiplex_region[0:cell_index_len]
                mt = multiplex_region[cell_index_len:]
                OUT.write("%s\t%s\t%s\n"%(read_id,cell_index,mt))
            else:
                j+=1

def extract_region(vector,error,cell_index_len,mt_len,isolator,read2_fastq,outfile,cores):
    '''
    A wrapper function to parallelize motif finding in the sequencing reads

    :param str vector: 5' flanking sequence for R2
    :param int error: number of snps/indels allowed in the vector sequence
    :param int cell_index_len: length of the cell index region
    :param int mt_len: length of the molecular tag region
    :param str read2_fastq: 2nd mate pair fastq location(can be .gz or not)
    :param str outfile: path to the output file
    :param int cores: number of CPUs to use
    :return: nothing
    :rtype:
    '''

    p = multiprocessing.Pool(cores)
    multiplex_len = cell_index_len + mt_len
    motif = r'((%s){e<=%i}([ACGT]{%i,%i})(ACG){s<=1}[ACGT]*)'%(vector,error,multiplex_len-1,multiplex_len+1)
    func = partial(find_motif,motif)
    print "\nLooking for the motif : %s in the sequencing reads.\n"%motif
    ## Print out results , takes in input a list of tuples which are processed in parallel by func()
    print_result(p.map(func,iterate_fastq(read2_fastq)),outfile,cell_index_len)
    p.close()
    p.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=
                                     argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-v','--vector',help="The vector nucleotide sequence",
                        required=True)
    parser.add_argument('--isolator',
                        help="The isolator sequence flanking the molecular tag",
                        required=True)
    parser.add_argument('-c',
                        '--cell-index-len',help="The length of the cell index",
                        type=int,required=True)
    parser.add_argument('-m','--mt-len',help="The length of the molecular tag",
                        type=int,required=True)
    parser.add_argument('-i','--read2-fastq',help="The input R2 fastq file",
                        required=True)
    parser.add_argument('-o','--out',help="The output",required=True)
    parser.add_argument('-e','--error',help="The max number of permissible "
                        "errors (substitions/indels) in the vector sequence",
                        required=False,type=int,default=3)
    parser.add_argument('-n','--parallel',
                         help="The number of cores to use",
                         default=1,type=int,required=False)
    args = parser.parse_args()
    print "\nRunning program with arguments : {0}\n".format(args)
    extract_region(args.vector,args.error,args.cell_index_len,args.mt_len,args.isolator,
                             args.read2_fastq,args.out,args.parallel)
