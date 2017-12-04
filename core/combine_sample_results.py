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

def clean_for_clustering(combined_cell_metrics_file,combined_umi_counts_file):
    ''' Clean the cells , i.e. remove cells with no ERCC reads

    :param: str combined_cell_metrics: the path to the combined metrics file
    :param: str combined_cell_metrics: the path to the combined umi counts file
    ''' 
    clean_cells = []
    with open(combined_cell_metrics_file,'r') as IN,open(combined_cell_metrics_file+'.clean','w') as OUT:
        for line in IN:
            line = line.strip('\n')
            if line.startswith('Cell'):
                print >> OUT,line
                continue
            else:                
                contents = line.split('\t')
                cell = contents[0]
                if int(contents[5]) == 0: ## Remove cells with no reads mapped to ERCC
                    continue
                else:
                    print >> OUT,line
                    clean_cells.append(cell)
                    
    with open(combined_umi_counts_file,'r') as IN,open(combined_umi_counts_file+'.clean','w') as OUT:
        for line in IN:
            line = line.strip('\n')
            contents = line.split('\t')
            if line.startswith('chromosome'):
                header_anno = contents[0:6]
                header_cells = contents[6:]
                outheader = '\t'.join(header_anno) + '\t' + '\t'.join(clean_cells)
                print >> OUT,outheader
                continue
            else:
                umis = contents[6:]
                umi_dict = dict(zip(header_cells,umis))
                ## Keep only the cells filtered in the metrics file
                out = contents[0:6]
                for cell in clean_cells:                     
                    out.append(umi_dict[cell])
                outline = '\t'.join(out)
                print >> OUT,outline
           
def sort_by_cell(outputfile):
    ''' Sort the output count files by Sample_Cells

    :param: str outputfile: path to the output file
    '''
    temp=outputfile+'.sorted'
    with open(outputfile,'r') as IN,open(temp,'w') as OUT:
        for line in IN:
            contents=line.strip('\n').split('\t')
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
        sample = os.path.basename(metric_file).strip("_read_stats.txt")
        assert sample != '', "Error could not identify sample name from file path : {}".format(metric_file)
        for line in IN:
            metric,val = line.strip('\n').split(':')
            metric_dict[metric][sample] = float(val)
    return metric_dict

def combine_sample_metrics(files_to_merge,outfile,is_low_input):
    ''' Combine metrics on the sample level similar to the cells
    :param list files_to_merge: the files to merge
    :param outfile: the output file to write to
    :param str: is_lowinput: Whether the protocol was for a low input application(1/0)
    '''
    sample_metrics = MyOrderedDict()
    ## Read metrics for each sample
    for sfile in files_to_merge:
        sample_metrics = read_sample_metrics(sfile,sample_metrics)
    ## Combine and write resultant output file
    with open(outfile,'w') as OUT:
        i=0
        for metric in sample_metrics:
            if metric=="num_reads_mapped_ercc" and is_low_input=="1":
                continue
            if i == 0:
                header = 'Samples\t'+'\t'.join(sample_metrics[metric].keys())
                OUT.write(header+'\n')
                i+=1
            out = metric
            for sample in sample_metrics[metric]:
                out = out+'\t'+float_to_string(round(sample_metrics[metric][sample],2))
            OUT.write(out+'\n')
        
def combine_cell_metrics(files_to_merge,outfile,cells_to_restrict):
    ''' Combine cell metrics from different samples
    :param list files_to_merge: the files to merge
    :param str outfile: the outputfile to write the aggregate metrics
    :param list cells_to_restrict: restrict cells to this list
    '''    
    cell_metrics = MyOrderedDict()
    files_to_merge  = natsort.natsorted(files_to_merge)
    for cfile in files_to_merge:
        cell_metrics = read_cell_file(cfile,cell_metrics)

    with open(outfile,'w') as OUT:
        i=0
        for cell in cell_metrics:
            if cell not in cells_to_restrict:
                continue
            if i == 0: ## Write Header
                header = 'Cells\t'+'\t'.join(cell_metrics[cell].keys())
                OUT.write(header+'\n')
                i+=1
            out = cell            
            for metric in cell_metrics[cell]:
                out = out+'\t'+cell_metrics[cell][metric]
            OUT.write(out+'\n')

def combine_count_files(files_to_merge,outfile,wts,cells_to_restrict=[]):
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
    :param list cells_to_restrict: Restrict cells to only this set when writing the primer count file
                                   This list is based on the criteria mentioned below
    
    :return The cells written to the file , any cell with < 5 UMIs for each gene is not written
    :rtype list
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
            cell_key = sample_name+'_'+str(cell)
            for line in IN:
                k1,k2,k3,k4,k5,k6,umi = line.rstrip('\n').split('\t')
                key = (k1,k2,k3,k4,k5,k6)
                ## Hash the umi counts by annotation and cells
                UMI[key][cell_key] = umi
                check_counts.append(int(umi))
                
            if not wts:
                if cell_key in cells_to_restrict:
                    header_cells.add(cell_key)
            else:
                if any(e >= 5 for e in check_counts): ## Check to make sure the cell has atleast 5 UMI count for any 1 gene
                    header_cells.add(cell_key)          
    ## Create header
    if wts:
        header = "chromosome\tstart\tstop\tstrand\tgene\tgene_type\t{cells}\n"
    else:
        header = "chromosome\tstart\tstop\tstrand\tgene\tprimer_sequence\t{cells}\n"
    temp = '\t'.join(list(header_cells))
    head = header.format(cells=temp)
    ## Print output
    with open(outfile,'w') as OUT:
        OUT.write(head)
        for key in UMI:
            write=True
            out = '\t'.join(key)
            for cell in header_cells:
                if cell not in UMI[key]:
                    raise Exception("Cell not hashed for Gene/Primer : {cell}-{k}".format(cell=cell,k=key))
                    #out = out + '\t0'
                else:
                    out = out + '\t{}'.format(UMI[key][cell])
            if write:        
                OUT.write(out+'\n')                
    ## Sort the count file
    sort_by_cell(outfile)

    return header_cells
