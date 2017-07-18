#cython: boundscheck=False, wraparound=False, nonecheck=False

cpdef bint compare_seq(str s1, str s2, int max_count):
    '''
    cython function to calculate to return whether two strings have a
    given max_mismatch
    '''
    cdef int i = 0
    cdef int count = 0
    cdef int l1
    cdef int l2
    l1 = len(s1)
    l2 = len(s2)
    
    if l1 < l2:
        bound = l1
    else:
        bound = l2

    while (i < bound and count <= max_count):
        if s1[i] != s2[i]:
            count+=1
        i+=1

    if count <= max_count:
        return True
    else:
        return False
