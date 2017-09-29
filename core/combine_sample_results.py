import sys
import glob
import os
import natsort
from collections import defaultdict,OrderedDict
from merge_mt_files import float_to_string

class MyOrderedDict(OrderedDict):
    def __missing__(self,key):
        val = self[key] = MyOrderedDict()
        return val

def sort__by_cell(outputfile):
    ''' Sort the output count files by Sample_Cells
    '''
    temp=outputfile+'.sorted'
    with open(outputfile,'r') as IN,open(temp,'w') as OUT:
        for line in IN:
            if contents[0] == 'chromosome': #header
                contents = line.strip('\n').split('\t')
                anno_header = '\t'.join(contents[0:6])
                cells = contents[6:]
                sorted_cells = '\t'.join(natsort.natsorted(cells))
                OUT.write(anno_header+'\t'+sorted_cells+'\n')
                continue
            anno = '\t'.join(contents[0:6])
            sorted_vals = [x for _,x in natsort.natsorted(zip(cells,contents[6:]))]
            OUT.write(anno+'\t'+'\t'.join(sorted_vals)+'\n')
    os.system('mv {temp} {outputfile}'.format(temp=temp,outputfile=outputfile))

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
            if contents[1:] == ['0']*len(contents[1:]): ## Skip cells with all zeros
                continue
            cell= contents[0]
            for i,metric in enumerate(metrics[1:]):
                metric_dict[cell][metric]= contents[i+1]

    return metric_dict

def read_sample_metrics(metric_file,metric_dict):
    ''' Read a sample metrics file
    '''
    with open(metric_file,'r') as IN:
        sample = os.path.basename(metric_file)
        assert sample != '', "Error could not identify sample name from file path : {}".format(metric_file)
        for line in IN:
            metric,val = line.strip('\n').split(':')
            metric_dict[sample][metric] = float(val)
    return metric_dict

def combine_sample_metrics(files_to_merge,outfile):
    ''' Combine metrics on the sample level similar to the cells
    :param list files_to_merge: the files to merge
    :param outfile: the output file to write to
    '''
    sample_metrics = MyOrderedDict()
    ## Read metrics for each sample
    for sfile in files_to_merge:
        sample_metrics = read_sample_metrics(sfile,sample_metrics)
    ## Combine and write resultant output file
    with open(outfile,'w') as OUT:
        i=0
        for sample in sample_metrics:
            if i == 0:
                header = 'Samples\t'+'\t'.join(sample_metrics[sample].keys())
                OUT.write(header+'\n')
                i+=1
            out = sample
            for metric in sample_metrics[sample]:
                out = out+'\t'+float_to_string(round(sample_metrics[sample][metric],2))
            OUT.write(out+'\n')
        
def combine_cell_metrics(files_to_merge,outfile):
    ''' Combine cell metrics from different samples
    :param list files_to_merge: the files to merge
    :param str outfile: the outputfile to write the aggregate metrics
    '''    
    cell_metrics = MyOrderedDict()
    for cfile in files_to_merge:
        cell_metrics = read_cell_file(cfile,cell_metrics)

    with open(outfile,'w') as OUT:
        i=0
        for cell in cell_metrics:
            if i == 0: ## Write Header
                header = 'Cells\t'+'\t'.join(cell_metrics[cell].keys())
                OUT.write(header+'\n')
                i+=1
            out = cell
            for metric in cell_metrics[cell]:
                out = out+'\t'+cell_metrics[cell][metric]
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

    :param list files_to_merge: full path to the files to merge
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
    head = header.format(cells=temp)
    print head
    ## Print output
    with open(outfile,'w') as OUT:
        OUT.write(head)
        for key in UMI:
            out = '\t'.join(key)
            for cell in header_cells:
                if cell not in UMI[key]:
                    out = out + '\t0'
                else:
                    out = out + '\t{}'.format(UMI[key][cell])
            OUT.write(out+'\n')      
                
