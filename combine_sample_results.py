import sys
import glob
import numpy as np
import pandas as pd
from collections import defaultdict,OrderedDict

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
                header = contents[6:]
                continue
            gene_info.append('\t'.join(contents[0:6]))
            counts.append(contents[6:])

    return (header,gene_info,counts)

def read_cell_file(cfile,metric_dict):
    '''
    '''

    i=0
    columns = []

    with open(cfile,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            if i==0: ## Header
                metrics = contents
                i+=1
                continue
            if contents[1:] == ['0']*len(contents[1:]):
                continue

            cell= contents[0]
            for i,metric in enumerate(metrics[1:]):
                metric_dict[metric][cell]= contents[i+1]

    return metric_dict

def combine_cell_metrics(output_folder,outfile):
    '''
    '''

    class MyOrderedDict(OrderedDict):
        def __missing__(self,key):
            val = self[key] = MyOrderedDict()
            return val
    cell_metrics = MyOrderedDict()

    cellfiles = glob.glob('{out}/*/*_cell_stats.txt'.format(out=output_folder))
    for cfile in cellfiles:
        cell_metrics = read_cell_file(cfile,cell_metrics)

    with open(outfile,'w') as OUT:
        i=0
        for metric in cell_metrics:
            print metric
            if i == 0: ## Write Header
                header = 'Metrics/Cells\t'+'\t'.join(cell_metrics[metric].keys())
                OUT.write(header+'\n')
                i+=1
            out = metric
            for cell in cell_metrics[metric]:
                out = out+'\t'+cell_metrics[metric][cell]
            OUT.write(out+'\n')

def combine_count_files(output_folder,outfile):
    ''' Combine multiple sample count files into 1
    '''

    prev_info = None
    prev_header = None
    samples = glob.glob('{out}/*/*.umi.counts.txt'.format(out=output_folder))
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
            count_matrix = np.empty(shape=(ngenes,nsamples*ncells),dtype='S9')
            print np.shape(count_matrix)
        headers.extend(header)

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
                Exception("Error in dimensions when assinging values to count matrix")
            i+=1

    with open(outfile,'w') as OUT:
        ## Convert numpy array to pandas dataframe for column name subsetting
        temp = pd.DataFrame(data=count_matrix,
                            index=gene_info,
                            columns=headers)
        print temp.shape
        temp = temp[temp.columns[(temp != '0').any()]]
        print temp.shape
        count_matrix = temp.as_matrix()
        new_headers = list(temp)

        annotation = "chrom\tstart\tstop\tstrand\tgene\tgene_type"
        head = '\t'.join(new_headers)

        OUT.write(annotation+'\t'+head+'\n')
        for i in range(0,len(gene_info)):
            out = gene_info[i] + '\t' + '\t'.join(count_matrix[i,:])+ '\n'
            OUT.write(out)


if __name__ == '__main__':
    #combine_count_files(sys.argv[1],sys.argv[2])
    combine_cell_metrics(sys.argv[1],sys.argv[3])
