import sys
import glob
import numpy as np

def read_count_file(countfile):
    ''' Read a count file
    '''

    gene_info = []
    counts = []
    header = []
    i = 0
    with open(countfile,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            if i == 0:
                i+=1
                header = '\t'.join(contents[6:])
                continue
            gene_info.append(contents[0:6])
            counts.append(contents[6:])

    return (header,gene_info,counts)

def combine_count_files(output_folder,outfile):
    ''' Combine multiple sample count files into 1
    '''

    prev_info = None
    prev_header = None
    samples = glob.glob('{out}/*/*.mt.counts.txt'.format(out=output_folder))
    nsamples = len(samples)
    prev_column = 0
    headers = []
    for countfile in samples:
        header,gene_info,counts = read_count_file(countfile)
        if prev_info:
            assert prev_info[0] == gene_info[0],"Mismatch in column ordering"
            assert prev_info[1] == gene_info[1],"Mismatch in column ordering"
            assert prev_info[-1] == gene_info[-1],"Mismatch in column ordering"
            prev_column = (prev_column + ncells)
        else:
            prev_info = gene_info
            prev_header = header
            ngenes = len(gene_info)
            ncells = len(counts[0])
            count_matrix = np.zeros(shape=(ngenes,nsamples*ncells))
            print np.shape(count_matrix)
        headers.append(header)
        ## Update count matrix
        i=0
        print countfile
        print prev_column,prev_column+ncells
        print ngenes
        bound = prev_column+ncells
        for c in counts:
            try:
                count_matrix[i,prev_column:bound] = c
            except:
                print c
                print count_matrix[i,:]
                Exception("zz")
            i+=1

    with open(outfile,'w') as OUT:
        annotation = "chrom\tstart\tstop\tstrand\tgene\tgene_type"
        head = '\t'.join(header)
        OUT.write(annotation+'\t'+head+'\n')
        for i in range(0,len(gene_info)):
            out = '\t'.join(gene_info[i]) + '\t' + '\t'.join([str(c) for c in count_matrix[i:]]) + '\n'
            OUT.write(out)


if __name__ == '__main__':
    combine_count_files(sys.argv[1],sys.argv[2])
