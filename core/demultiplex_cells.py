import sys
import gzip
import io
import os
from extract_multiplex_region import open_by_magic

"""
Demultiplex the R1 fastq file into individual cells based
upon the cell index, incorporate the barcode tag to the ReadID
for downstream anaylsis.
"""


def iterate_fastq():
    '''
    '''


def read_cell_index_file(cell_index_file):
    '''
    Read the file containing the cell indices and
    return it as a list
    :param str cell_index_file: the cell index file
    :return: list containing the cell indices
    :rtype: list of str
    '''
    with open(cell_index_file,'r') as IN:
        return [line.strip('\n') for line in IN]

def create_cell_dirs(base_dir,cell_indices):
    ''' Create directories to store cell specific read info
    :param list cell_indices: list containing the cell indices
    '''
    map(lambda x:os.makedirs(os.path.join(
        base_dir,x[1]+'_cell_'+str(x[0]))),
        enumerate(cell_indices))

cell_indices = read_cell_index_file(sys.argv[1])
create_cell_dirs('./',cell_indices)
