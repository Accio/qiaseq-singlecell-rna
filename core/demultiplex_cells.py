import sys
import gzip
import io
import os
import errno
import editdistance
import collections
import regex
import re
from extract_multiplex_region import open_by_magic

"""
Demultiplex the R1 fastq file into individual cells based
upon the cell index, incorporate the barcode tag to the ReadID
for downstream anaylsis.
"""

## To do :
## 1. Improve runtime, currently takes a ~64 secs to multiplex and write fastqs
## for a readset of size 2 million reads


def benchmark(func):
    '''
    A decorator that prints the time a function takes
    to execute.
    :param func func: function to decorate
    :return: a function object
    :rtype: func
    '''
    import time
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        print func.__name__, time.clock()-t
        return res

    return wrapper

def trim_read(read,seq_pattern):
    ''' Trim a read for a given sequence, will make this function more generic
    in the future , for now this will simply chop off the polyA tail from R1
    reads , for e.g. AAGTCTGGCCATGCAAAAAAAAA ---> AAGTCTGGCCAT

    :param str read: the read sequence
    :param re object seq_pattern: regular expression to match the polyA tail
    :return: a tuple of the trimmed sequence and the start position of the trimmed region
    :rtype: tuple
    '''

    match = seq_pattern.match(read)
    if match:
        trimmed_seq = match.group(1)
        tail = match.group(2)
        #print "Trimming {polyT} from the read {read}".format(polyT=tail,read=read)
        return (trimmed_seq,match.start(2))
    else:
        return (read,0)

def iterate_fastq(fastq):
    '''
    Read a fastq file and return the 4 lines as a list
    '''

    with open_by_magic(fastq) as IN:
        while True:
            yield [IN.next().rstrip('\n'),IN.next().rstrip('\n'),
                   IN.next().rstrip('\n'),IN.next().rstrip('\n')]

def read_multiplex_file(multiplex_file):
    '''
    Read the mutliplex file created using extract_multiplex_region.py

    :param str multiplex_file: path to the multiplex region file
    :yield: list of (read2_id,cell_index,mt)
    '''

    with open(multiplex_file,'r') as IN:
        IN.readline()
        for line in IN:
            yield line.rstrip('\n').split('\t')

def match_cell_index(cell_indices,cell_indices_mismatch,cell_index,edit_dist):
    '''
    Edit distance match on the cell index

    :param dict cell_indices: {cell_index:cellnumber}
    :param dict cell_indices_mismatch: {mismatched_cell_index:cell_index}
    :param str cell_index
    :param int edit_dist
    :return: Whether the cell index matches an valid list of indices
    :rtype: Bool
    '''
    if cell_index in cell_indices:
        return (True,cell_index)
    else:
        if edit_dist == 1:
            if cell_index in cell_indices_mismatch:
                return (True,cell_indices_mismatch[cell_index])
            else:
                return (False,None)
        else:
            raise Exception("Only 1 mismatch allowed in cell index !")

def create_cell_index_db(multiplex_file):
    '''
    Create a db to store read_id -> cell_index, mt
    Use this for out of memory situations , use the
    shelve module in python.
    '''

def create_read_id_hash(multiplex_file):
    '''
    Creates a dict read_id -> [cell_index,mt]
    :param str multiplex_file: a tsv <read_id> <cell_index> <mt>
    :return: a dictionary of {read_id:[cell_index,mt]}
    :rtype: dict
    '''
    d = {}
    for read_id,cell_index,mt in read_multiplex_file(multiplex_file):
        key = read_id.split()[0]
        d[key] = [cell_index,mt]
    return d

def mutate(x):
    ''' Returns all possible single base substitutions of a dna string
    including N , to represent sequencing error
    :param str x: the dna nucleotide sequence to mutate
    :yields the possible mutated strings
    '''
    substitions = ['A','C','G','T','N']
    for i in xrange(len(x)):
        for b in substitions:
            temp = list(x)
            if temp[i] != b:
                temp[i] = b
                yield "".join(temp)

def read_cell_index_file(cell_index_file):
    '''
    Read the file containing the cell indices and
    return it as a dict

    :param str cell_index_file: the cell index file
    :return: tuple of dicts
    :rtype: (dict,dict)
    :raises: Exception for duplicate cell index
    '''                
    d = collections.OrderedDict()
    mismatch_d = {}
    i=1
    with open(cell_index_file,'r') as IN:
        for line in IN:
            key = line.rstrip('\n')
            if key in d:
                raise Exception('Duplicate cell index encountered !')
            d[key] = i
            i+=1
            # all possible single base mutations
            for m_key in mutate(key):
                if m_key in mismatch_d:
                    raise Exception("Duplicate mutated cell index encountered !")
                mismatch_d[m_key] = key
                
    return (d,mismatch_d)

def write_metrics(metric_file,metric_dict,metrics):
    ''' Write Metrics
    :param str metric_file: output file to write the metrics to
    :param dict metric_dict: dictionary containing the metrics
    :param list metrics: ordered list of metrics containing dictionary keys
    '''
    with open(metric_file,'w') as OUT:
        for key in metrics:
            if key in metric_dict:
                val = metric_dict[key]
            else:
                val = 0        
            OUT.write('{metric}: {value}\n'.format(metric=key,value=val))


def write_fastq(read_info,OUT):
    ''' Write a fastq file

    :param str read_info: a list whose elements are the 4 lines of a fastq
    :param fastq_loc: full path to the fastq file to write to
    :return: nothing
    '''

    #out = '\n'.join(read_info)
    OUT.write(read_info+'\n')

@benchmark
def create_cell_fastqs(base_dir,metric_file,cell_index_file,
                       cell_multiplex_file,read_file1,editdist,wts=False):
    '''
    Demultiplex and create individual cell fastqs

    :param str base_dir: base directory to create subfolders and files
    :param str metric_file: name of the metric file to write to
    :param str cell_index_file: file containing the valid cell indices
    :param str cell_multiplex_file: tsv file <read2_id> <cell_index> <mt>
    :param str read_file1: read1 fastq file
    :param int editdist: whether to allow a single base mismatch in the cellindex. 0 or 1
    :param bool wts: whether this is for wts
    :return: nothing
    '''

    ## Initialize variables
    FASTQS= {}
    METRICS = {}
    num_reads=0
    reads_demultiplexed=0
    reads_dropped_cellindex=0
    reads_dropped_size = 0
    if wts:
        polyA_motif = re.compile(r'^([ACGTN]*?[CGTN])([A]{9,}[ACGNT]*$)')
    else:
        polyA_motif = re.compile(r'^([ACGTN]{42,}[CGTN])([A]{8,}[ACGNT]{1,}$)')
    cell_indices,cell_indices_mismatch = read_cell_index_file(cell_index_file)
    read_id_hash = create_read_id_hash(cell_multiplex_file)
    reads_to_demultiplex = len(read_id_hash.keys())
    ## Create cell specific dirs and open file handles
    def mkdir_p(path):
        try:
            os.makedirs(path)
        except OSError as exc: ## Hnadle race condition
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise exc
    map(lambda x:mkdir_p(os.path.join(base_dir,'Cell'+str(x[0]+1)+'_'+x[1])),
        enumerate(cell_indices))

    for cell_index,cell_num in cell_indices.items():
        fastq=os.path.join(base_dir,'Cell'+str(cell_num)+'_'+cell_index+
                                         '/cell_'+str(cell_num)+'_R1.fastq')
        metric=os.path.join(base_dir,'Cell'+str(cell_num)+'_'+cell_index+
                                         '/cell_'+str(cell_num)+'_demultiplex_stats.txt')
        FASTQS[cell_index] = open(fastq,'w')
        METRICS[cell_index] = metric
    ## Store metrics for each cell
    cell_metrics = collections.defaultdict(lambda:collections.defaultdict(int))

    ## Iterate over the R1 fastq , check if the cell index is valid,
    ## 3' polyA tail and write as a fastq file
    for read_id,seq,p,qual in iterate_fastq(read_file1):
        key = read_id.split()[0]
        if key in read_id_hash:
            cell_index,mt = read_id_hash[key]
            match,cell_index = match_cell_index(cell_indices,cell_indices_mismatch,cell_index,editdist)
            if match == True:
                new_read_id = key+":%s"%mt ## CellIndex can be added here as well
                trimmed_seq,trimming_index = trim_read(seq,polyA_motif)
                if trimming_index: ## Non zero trimming index
                    qual = qual[0:trimming_index]
                if len(trimmed_seq) < 25: ## Drop reads below 25 b.p
                    reads_dropped_size+=1
                    num_reads+=1
                    cell_metrics[cell_index]['reads total']+=1
                    continue
                read_info = new_read_id+'\n'+trimmed_seq+'\n'+p+'\n'+qual
                write_fastq(read_info,FASTQS[cell_index])
                cell_metrics[cell_index]['reads total']+=1
                cell_metrics[cell_index]['after_qc_reads']+=1
                reads_demultiplexed+=1
            else:
                reads_dropped_cellindex+=1
        num_reads+=1

    for cell in FASTQS:
        FASTQS[cell].close()

    ## Write out the metrics
    ## 1. Per cell level
    metrics_to_write = ['reads total','after_qc_reads']
    for cell,mfile in METRICS.items():
        write_metrics(mfile,cell_metrics[cell],metrics_to_write)
    ## 2. Aggregated across all the cells
    metric_dict = collections.OrderedDict(
        [('reads total',num_reads),
         ('reads dropped, cell id not extracted',num_reads - reads_to_demultiplex),
         ('reads dropped, cell id not matching a used oligo within edit distance {} bp'.format(editdist),reads_dropped_cellindex),
         ('reads dropped, less than 25 bp',reads_dropped_size)
        ]
    )
    write_metrics(metric_file,metric_dict,metric_dict.keys())

if __name__ == '__main__':
    create_cell_fastqs(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])


