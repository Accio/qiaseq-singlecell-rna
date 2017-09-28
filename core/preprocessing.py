import re
import os
from extract_multiplex_region import open_by_magic

def random_sample_fastq():
    ''' Using Algorithm R(Reservoir Sampling)
    https://en.wikipedia.org/wiki/Reservoir_sampling
    O(n) with O(1) memory foot print
    Can't do better than O(n) if we want a true random
    sampling of the data. This algorithm will scale well
    with size of the fastq.
    '''    
    ## The random line chosen will be the read_id line of a fastq
    rand_lines = {} ## Random lines hashed by read_id
    

def is_vector_present(fastq,vector):
    ''' Function to identify whether a vector is present or not
    :param str fastq: the input R1/R2 fastq file
    :param str vector: the vector sequence to identify
    :param str instrument: instrument used for sequencing, used for determining 
                           how to identify the vector, i.e. quality or sequence
    '''
    pass

def preprocess_reads(R1,R2,vector):
    '''
    '''
    ## Initialize some variables
    temp_R1 = R1+'.bak'
    temp_R2 = R2+'.bak'

    ## Sample randomly from R2 reads
    
    if is_vector_present(R1,vector,instrument):
        trim_vector(R1,temp_R1,vector)
    if not is_vector_present(R2,vector,instrument):
        add_vector(R2,temp_R2,vector)

    ## Rename files 
    os.system('mv {temp_R1} {R1}'.format(temp_R1=temp_R1,R1=R1))
    os.system('mv {temp_R2} {R2}'.format(temp_R2=temp_R2,R2=R2))
    return instrument
