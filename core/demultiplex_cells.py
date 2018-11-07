from pathos import multiprocessing
import gzip
import functools
import logging
import sys
import collections
import os
import errno
import re
import edlib

import regex

import pyximport
pyximport.install(reload_support=True)
from _utils import two_fastq_heads

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
logger.addHandler(ch)


# Metric names
# 1. Per cell level
CELL_READS_TOTAL = "reads total"
CELL_AFTER_QC    = "after_qc_reads" # for now any read >= 25 b.p after polyA trim
# 2. Sample Index level
OVERALL_TOTAL                         =  "reads total"
OVERALL_DROPPED_ALL_N                 =  "reads dropped, all NNNNNN sequence"
OVERALL_DROPPED_CELLID_NOT_EXTRACTED  =  "reads dropped, cell id not extracted"
OVERALL_DROPPED_CELLID_MISMATCH       =  "reads dropped, cell id not matching a used oligo within edit distance {e} bp"
OVERALL_DROPPED_LT_25BP               =  "reads dropped, less than 25 bp"


def open_fh(fname1,fname2,read=True):
    ''' Return appropriate file handles
    :param str fname1: R1 fastq file name
    :param str fname2: R2 fastq file name
    :param bool read: rb/wb mode
    :rtype tuple
    :returns tuple of file handles
    '''
    mode = "rb" if read else "wb"
    if fname1.endswith(".gz"):
        return (gzip.open(fname1,mode),gzip.open(fname2,mode))
    else:
        return (open(fname1,mode),open(fname2,mode))

def close_fh(fh1,fh2):
    ''' Closes file handles
    :param str fh1: R1 file handle
    :param str fh2: R2 file handle
    '''
    fh1.close()
    fh2.close()
    
def iterate_fastq(f,f2,ncpu,buffer_size=4*4*1024**2):
    ''' Copied from cutadapt, 
    added logic to yield a list of buffers equal to the number of CPUs
    :param file_handle f:  R1 fastq
    :param file_handle f2: R2 fastq
    :param int ncpu: length of buffer list 
    :param int buffer_size: size in bytes of each element of buffer list

    :yields list of length ncpu
    '''
    buf1 = bytearray(buffer_size)
    buf2 = bytearray(buffer_size)
    
    # Read one byte to make sure we are processing FASTQ
    start1 = f.readinto(memoryview(buf1)[0:1])
    start2 = f2.readinto(memoryview(buf2)[0:1])

    if (start1 == 1 and buf1[0:1] != b'@') or (start2 == 1 and buf2[0:1] != b'@'):
        raise Exception('Paired-end data must be in FASTQ format when using multiple cores')

    to_yield = []
    nchunks = 0
    while True:
        bufend1 = f.readinto(memoryview(buf1)[start1:]) + start1
        bufend2 = f2.readinto(memoryview(buf2)[start2:]) + start2
        if start1 == bufend1 and start2 == bufend2:
            break

        end1, end2 = two_fastq_heads(buf1, buf2, bufend1, bufend2)
        assert end1 <= bufend1
        assert end2 <= bufend2

        if end1 > 0 or end2 > 0:
            nchunks+=1
            to_yield.append((memoryview(buf1)[0:end1].tobytes(), memoryview(buf2)[0:end2].tobytes()))
            
        start1 = bufend1 - end1
        assert start1 >= 0
        buf1[0:start1] = buf1[end1:bufend1]
        start2 = bufend2 - end2
        assert start2 >= 0
        buf2[0:start2] = buf2[end2:bufend2]
        
        if nchunks == ncpu:
            yield to_yield
            to_yield = []
            nchunks = 0

    if start1 > 0 or start2 > 0:
        to_yield.append((memoryview(buf1)[0:start1].tobytes(), memoryview(buf2)[0:start2].tobytes()))
    if len(to_yield) > 0:
        yield to_yield

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

def read_cell_index_file(cell_index_file,cell_indices_used):
    '''
    Read the file containing the cell indices and
    return cell indices used and cell indices different by 1 b.p

    :param str cell_index_file: the cell index file
    :param str cell_indices_used: comma delimeted cell ids used in the experiment
    :return: tuple of dicts
    :rtype: (dict,dict)
    :raises: Exception for duplicate cell index
    '''                
    d = collections.OrderedDict()
    mismatch_d = {}
    i=1
    used_indices = set(cell_indices_used.split(','))
    with open(cell_index_file,'r') as IN:
        for line in IN:
            key = line.rstrip('\n')
            if key in d:
                raise Exception('Duplicate cell index encountered !')
            if 'C'+str(i) in used_indices or cell_indices_used == 'all':            
                d[key] = i            
                # all possible single base mutations
                for m_key in mutate(key):
                    if m_key in mismatch_d:
                        raise Exception("Duplicate mutated cell index encountered !")
                    mismatch_d[m_key] = key
            i+=1                
    return (d,mismatch_d)
        
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: ## Hnadle race condition
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise exc

def id_cell_umi_nextseq(r2_seq,cell_index_len,umi_len):
    ''' Identify cell id and umi region for nextseq reads
    '''
    if len(r2_seq) >= cell_index_len + umi_len - 1: # allow 1 base offset
        cellid = r2_seq[0:cell_index_len]
        umi    = r2_seq[cell_index_len:cell_index_len+umi_len]
        return (cellid,umi)
    else:
        return (None,None)

def id_cell_umi(r2_seq,cell_index_len,umi_len,vector,error):
    ''' Identify cell id and umi region for MiSeq/HiSeq reads
    '''
    found = False
    alignment = edlib.align(vector,r2_seq,mode="SHW")    
    if alignment["editDistance"] <= error:   
    
        vector_end_pos = alignment["locations"][-1][1] # 0-based position on r2
        temp = r2_seq[vector_end_pos + 1:]

        if len(temp) >= cell_index_len + umi_len - 1: # allow 1 base offset
            found = True

            cellid = temp[0:cell_index_len]
            umi    = temp[cell_index_len:cell_index_len+umi_len]       

    if found:
        return (cellid,umi)
    else:
        return (None,None)
    
    
def process_reads(args,buffer_):
    ''' Process R1,R2 fastq files , identify and trim synthetic oligos and polyA tail
    '''
    cell_indices,cell_indices_mismatch,editdist,wts,cell_index_len,umi_len,vector,error,instrument = args
    
    # unpack input byte string                                 
    buff_r1,buff_r2 = buffer_
    r1_lines        = buff_r1.split(b"\n")
    r2_lines        = buff_r2.split(b"\n")

    # init counters
    num_reads                               = 0
    reads_dropped_all_N                     = 0
    reads_dropped_cellid_not_extracted      = 0
    reads_dropped_cellid_not_matching_oligo = 0
    reads_dropped_lt_25bp                   = 0

    out_lines_r1 = collections.defaultdict(list)
    cell_metrics = collections.defaultdict(
        lambda:collections.defaultdict(int))
    
    i = 1
    for line in zip(r1_lines,r2_lines):
        
        if line[0] == "": # last element is empty because of the split("\n") above
            continue
        if i % 4 == 1: # header
            r1_readid,r2_readid = line
        elif i % 4 == 2: # seq
            r1_seq,r2_seq = line
        elif i % 4 == 3: # placeholder
            pass
        elif i % 4 == 0: # qual
            num_reads += 1
            r1_qual,r2_qual = line

            # have R1 and R2 ready to process now
            temp_r1_readid = r1_readid.split(" ")[0]
            temp_r2_readid = r2_readid.split(" ")[0]
            if temp_r1_readid != temp_r2_readid:
                raise UserWarning("demultiplex_cells:R1,R2 read ids are not in sync : R1:{r1} ; R2:{r2}".format(r1=temp_r1_readid,r2=temp_r2_readid))
            elif len(r1_qual) != len(r1_seq) or len(r2_qual) != len(r2_seq):
                raise UserWarning("demultiplex_cells:Read has different length qual and seq strings {}}".format(temp_r1_readid))

            # extract cell id , umi region from R2
            stop = False
            if instrument.upper() == "NEXTSEQ":
                cellid,umi = id_cell_umi_nextseq(r2_seq,cell_index_len,umi_len)
            else:
                cellid,umi = id_cell_umi(r2_seq,cell_index_len,umi_len,vector,error)
                
            if cellid:
                if cellid not in cell_indices:                    
                    stop = True
                    if editdist == 1:
                        if cellid in cell_indices_mismatch: # if within 1 edit distance recover                        
                            cellid = cell_indices_mismatch[cellid]
                            stop = False
                        else:
                            reads_dropped_cellid_not_matching_oligo += 1
                    else:
                        reads_dropped_cellid_not_matching_oligo += 1 
            else:
                stop = True                
                if len(set(r2_seq)) == 1 and r2_seq[0] == 'N':
                    reads_dropped_all_N += 1
                else:
                    reads_dropped_cellid_not_extracted += 1

            if stop:
                i+=1                
                continue

            new_read_id = temp_r1_readid + b":{}".format(umi)

            # polyA trim
            polyA_pattern = _REGEX_["r1_polyA"]
            match = polyA_pattern.match(r1_seq)
            if match:
                tail = match.group(2)
                trimming_index  = match.start(2)                
                r1_trimmed_seq  = match.group(1)
                r1_trimmed_qual = r1_qual[0:trimming_index]
            else:
                r1_trimmed_seq  = r1_seq
                r1_trimmed_qual = r1_qual

            # store per cell level metrics
            cell_metrics[cellid][CELL_READS_TOTAL]+=1
            if len(r1_trimmed_seq) < 25:
                reads_dropped_lt_25bp+=1
                stop = True
            else:
                cell_metrics[cellid][CELL_AFTER_QC]+=1
                
            if stop:
                i+=1
                continue

            trimmed_r1_lines = b"\n".join([new_read_id,r1_trimmed_seq,b"+",r1_trimmed_qual])
            out_lines_r1[cellid].append(trimmed_r1_lines)
            
        i+=1

    metrics = (cell_metrics, num_reads, reads_dropped_all_N, reads_dropped_cellid_not_extracted, reads_dropped_cellid_not_matching_oligo, reads_dropped_lt_25bp)
    return (out_lines_r1,metrics)

def compile_regex(wts,instrument,multiplex_len,vector,error):
    ''' Setup regular expressions for cell id umi identification and polyA trim
    :param bool wts : Whether this is a polyA wts experiment
    :param str instrument : MiSeq/HiSeq or NextSeq
    :param int multiplex_len : CellID + UMI length
    :param str vector : Vector sequence
    :param int error : Indel/SNPs to tolerate on vector sequence    
    '''
    # regex used for trimming
    if wts:
        r1_polyA = re.compile(r'^([ACGTN]*?[CGTN])([A]{9,}[ACGNT]*$)')
    else:
        r1_polyA = re.compile(r'^([ACGTN]{42,}[CGTN])([A]{8,}[ACGNT]{1,}$)')

    # deprecated not used anymore; kept here for reference
    if instrument.upper() == 'NEXTSEQ':
        r2_structure = regex.compile(r'(([ACGTN]{%i})[ACGTN]*)'%(multiplex_len))
        match_group = 2
    else:
        match_group = 3
        r2_structure = regex.compile(r'((%s){e<=%i}([ACGNT]{%i,%i})[ACGTN]*)'%(
            vector,error,multiplex_len-1,multiplex_len+1))

    global _REGEX_
    _REGEX_ = {"r1_polyA":r1_polyA,"r2_structure":(r2_structure,match_group)}

def write_metrics(metric_file,metric_dict,metrics):
    ''' Write Metrics
    :param str metric_file: output file to write the metrics to
    :param dict metric_dict: dictionary containing the metrics
    :param list metrics: ordered list of metrics containing dictionary keys
    '''
    with open(metric_file,"w") as OUT:
        for key in metrics:
            if key in metric_dict:
                val = metric_dict[key]
            else:
                val = 0
            OUT.write("{metrict}: {value}\n".format(metrict=key, value=val))
    
def demux(r1,r2,cell_index_file,base_dir,out_metric_file,cell_indices_used,vector,
          instrument,wts,return_demux_rate,cell_index_len,umi_len,editdist,error,ncpu,buffer_size,
          verbose = False):
    ''' Demultiplex and write fastq files for each cell
    :param str r1: R1 fastq file
    :param str r2: R2 fastq file
    :param str cell_index_file: File with 1 line for each cell index oligo
    :param str base_dir: Base output directory
    :param str out_metric_file: Output metric file for overall sample index level metrics
    :param str cell_indices_used: Cell indices used in the experiment , i.e. C1,C2,C3,C4, etc.
    :param str vector : 25-mer on R1 5' end for MiSeq/HiSeq reads
    :param str instrument : MiSeq/HiSeq or NextSeq
    :param bool wts : Whether this is a polyA wts experiment
    :param bool return_demux_rate : Return reads_after_demux/total_reads if requested
    :param int cell_index_len : Cell Index oligo length
    :param int umi_len : UMI oligo length
    :param int editdist : 0/1 ; Whether to allow 1 editdistance mismatch to cell indices specified
    :param int error : Number of Indels/SNPs to tolerate in vector sequence
    :param int ncpu: Number of CPUs to use
    :param int buffer_size : Read these many MegaBytes(MB) from fastq file to memory for each read pair for each cpu
    :param bool verbose : Whether to verbosely log
    '''

    wts                = bool(wts)
    return_demux_rate  = bool(return_demux_rate)
    cell_index_len     = int(cell_index_len)
    umi_len            = int(umi_len)
    editdist           = int(editdist)
    error              = int(error)
    ncpu               = int(ncpu)
    buffer_size        = int(buffer_size)*1024**2

    FASTQS = {}
    METRICS = {}
    
    assert instrument.upper() in ["NEXTSEQ","MISEQ/HISEQ"], "Incorrect instrument specification"
    
    multiplex_len = cell_index_len + umi_len
    compile_regex(wts,instrument,multiplex_len,vector,error)
    cell_indices,cell_indices_mismatch = read_cell_index_file(cell_index_file,cell_indices_used)
   
    for cell_index,cell_num in cell_indices.items():
        path = os.path.join(base_dir,'Cell'+str(cell_num)+'_'+cell_index)
        mkdir_p(path)
        
    for cell_index,cell_num in cell_indices.items():
        fastq=os.path.join(base_dir,'Cell'+str(cell_num)+'_'+cell_index+
                                         '/cell_'+str(cell_num)+'_R1.fastq')
        metric=os.path.join(base_dir,'Cell'+str(cell_num)+'_'+cell_index+
                                         '/cell_'+str(cell_num)+'_demultiplex_stats.txt')
        FASTQS[cell_index] = open(fastq,'w')
        METRICS[cell_index] = metric

    
    f,f2 = open_fh(r1,r2)
    
    p = multiprocessing.Pool(ncpu)

    args = (cell_indices,cell_indices_mismatch,editdist,wts,cell_index_len,umi_len,vector,error,instrument)
    func = functools.partial(process_reads,args)
    
    nchunk                                   =  0
    total_reads                              =  0
    reads_dropped_all_N                      =  0
    reads_dropped_cellid_not_extracted       =  0
    reads_dropped_cellid_not_matching_oligo  =  0
    reads_dropped_lt_25bp                    =  0

    cell_metrics = collections.defaultdict(
        lambda:collections.defaultdict(int))    
    
    for chunks in iterate_fastq(f,f2,ncpu,buffer_size):
        res = p.map(func,chunks)
        for trimmed_r1_lines,metrics in res:
            
            # unpack return variables and update counters
            temp_cell_metrics = metrics[0]
            total_reads                             += metrics[1]
            reads_dropped_all_N                     += metrics[2]
            reads_dropped_cellid_not_extracted      += metrics[3]
            reads_dropped_cellid_not_matching_oligo += metrics[4]
            reads_dropped_lt_25bp                   += metrics[5]
            
            for cell in temp_cell_metrics: # accumulate cell specific                
                for metric in temp_cell_metrics[cell]:
                    cell_metrics[cell][metric] += temp_cell_metrics[cell][metric]

                out_r1 = b"\n".join(trimmed_r1_lines[cell])
                FASTQS[cell_index].write(out_r1)
                FASTQS[cell_index].write(b"\n")
                
        nchunk += 1
        if verbose:
            logger.info("Processed {} read fragments".format(total_reads))

    # write final metrics
    # 1. Per cell level
    metrics_to_write = [CELL_READS_TOTAL, CELL_AFTER_QC]
    for cell,mfile in METRICS.items():
        write_metrics(mfile,cell_metrics[cell],metrics_to_write)

    # 2. On a Sample Index level
    global OVERALL_DROPPED_CELLID_MISMATCH
    OVERALL_DROPPED_CELLID_MISMATCH = OVERALL_DROPPED_CELLID_MISMATCH.format(e = editdist)
    metric_dict = collections.OrderedDict([(OVERALL_TOTAL, total_reads),
                                           (OVERALL_DROPPED_ALL_N, reads_dropped_all_N),
                                           (OVERALL_DROPPED_CELLID_NOT_EXTRACTED, reads_dropped_cellid_not_extracted),
                                           (OVERALL_DROPPED_CELLID_MISMATCH, reads_dropped_cellid_not_matching_oligo),
                                           (OVERALL_DROPPED_LT_25BP, reads_dropped_lt_25bp)])
    write_metrics(out_metric_file, metric_dict, metric_dict.keys())

    # close file handles
    for cell in METRICS:
        FASTQS[cell].close()
    close_fh(f,f2)

    if return_demux_rate:
        passed_reads = float(total_reads - reads_dropped_all_N - \
                             reads_dropped_cellid_not_extracted - \
                             reads_dropped_cellid_not_matching_oligo - \
                             reads_dropped_lt_25bp)
        return passed_reads/total_reads
        
if __name__ == '__main__':
    demux(*sys.argv[1:])
