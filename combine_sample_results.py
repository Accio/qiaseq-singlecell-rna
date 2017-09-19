import sys
import glob
import os
from collections import defaultdict,OrderedDict

def read_cell_file(cfile,metric_dict):
    ''' Read a cell metrics file and return the parsed metrics as a dict

    :param str cfile: the cell metric file
    :param dict metric_dict: a dict of dict of metrics
    :return the dictionary of parsed metrics
    :rtype: dict
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

def combine_cell_metrics(files_to_merge,outfile):
    ''' Combine cell metrics from different samples
    :param str files_to_merge: the files to merge
    :param str outfile: the outputfile to write the aggregate metrics
    '''
    class MyOrderedDict(OrderedDict):
        def __missing__(self,key):
            val = self[key] = MyOrderedDict()
            return val
    cell_metrics = MyOrderedDict()

    for cfile in files_to_merge:
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

def combine_count_files(files_to_merge,outfile,wts):
    ''' Function to combine cells from different samples into 1 file
    The directory strucuture of a typical run loooks like : 
              <output_folder> 
                --- Sample1
                  --- cell1
                    --- cell1_Sample1
                --- Sample2
                    --- cell1
                      --- cell1_Sample2

    :param str files_to_merge: full path to the files to merge
    :param str outfile: The combined outputfile to write to
    :param bool wts: Whether this was whole transcriptome sequencing 
    '''
    i = 0
    UMI = defaultdict(lambda:defaultdict(int))
    header_cells = set()
    ## Iterate over the files to merge
    for f in files_to_merge:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        sample_name = os.path.dirname(f).split('/')[-2]
        check_counts = []
        with open(f,'r') as IN:
            cell_key = sample_name+'_Cell'+str(cell)
            for line in IN:
                k1,k2,k3,k4,k5,k6,umi = line.rstrip('\n').split('\t')
                key = (k1,k2,k3,k4,k5,k6)
                ## Hash the umi counts by annotation and cells
                UMI[key][cell_key] = umi
                check_counts.append(int(umi))
            ## Add check here to see if cells are part of user's listing
            if not all(e == 0 for e in check_counts): ## Check to make sure the cell doesnt have zero counts
                header_cells.add(cell_key)
    ## Create header
    if wts:
        header = "chromosome\tstart\tstop\tstrand\tgene\tgene_type\t{cells}\n"
    else:
        header = "chromosome\tstart\tstop\tstrand\tgene\tprimer_sequence\t{cells}\n"
    temp = '\t'.join(list(header_cells))
    header.format(cells=temp)
    ## Print output
    with open(outfile,'w') as OUT:
        OUT.write(header)
        for key in UMI:
            out = '\t'.join(key)
            for cell in header_cells:
                if cell not in UMI[key]:
                    out = out + '\t0'
                else:
                    out = out + '\t{}'.format(UMI[key][cell])
            OUT.write(out+'\n')
        
                
if __name__ == '__main__':
    combine_count_files(sys.argv[1],sys.argv[2])
    combine_cell_metrics(sys.argv[1],sys.argv[3])
