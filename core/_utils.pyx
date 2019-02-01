cimport cython
ctypedef fused bytes_or_bytearray:
    bytes
    bytearray

'''
These are cython functions I took from cutadapt's source code :
https://github.com/marcelm/cutadapt/tree/master/src/cutadapt
All credit for these functions are atributed to its author
'''

def two_fastq_heads(bytes_or_bytearray buf1, bytes_or_bytearray buf2, Py_ssize_t end1, Py_ssize_t end2):
    '''
    Skip forward in the two buffers by multiples of four lines.

    Return a tuple (length1, length2) such that buf1[:length1] and
    buf2[:length2] contain the same number of lines (where the
    line number is divisible by four).
    '''
    cdef:
        Py_ssize_t pos1 = 0, pos2 = 0
        Py_ssize_t linebreaks = 0
        unsigned char * data1 = buf1
        unsigned char * data2 = buf2
        Py_ssize_t record_start1 = 0
        Py_ssize_t record_start2 = 0

    while True:
        while pos1 < end1 and data1[pos1] != '\n':
            pos1 += 1
        if pos1 == end1:
            break
        pos1 += 1
        while pos2 < end2 and data2[pos2] != '\n':
            pos2 += 1
        if pos2 == end2:
            break
        pos2 += 1
        linebreaks += 1
        if linebreaks == 4:
            linebreaks = 0
            record_start1 = pos1
            record_start2 = pos2

    # Hit the end of the data block
    return record_start1, record_start2


def quality_trim(str qualities, str bases, int cutoff_back, bint is_nextseq, int base=33):
    '''
    '''
    cdef int s
    cdef int max_qual
    cdef int stop = len(qualities)
    cdef int start = 0
    cdef int i

    # 3' end trimming
    max_qual = 0
    s = 0
    for i in reversed(xrange(len(qualities))):
        if is_nextseq and bases[i] == "G":
            q = cutoff_back - 1
        else:
            q = ord(qualities[i]) - base
        s += cutoff_back - q

        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            stop = i

    return stop
