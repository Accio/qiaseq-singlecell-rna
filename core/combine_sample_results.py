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
    ''' Clean the header line in the output file for clustering analysis

    :param: str combined_cell_metrics: the path to the combined metrics file
    :param: str combined_cell_metrics: the path to the combined umi counts file
    ''' 
    clean_cells = []
    clean_header_metrics = ["reads_total","reads_used_aligned_to_genome","reads_used_aligned_to_ERCC","UMIs","detected_genes"]
    clean_header_umi = ["gene_id","gene","strand","chrom","loc_5prime_grch38","loc_3prime_grch38"]        
    
    with open(combined_cell_metrics_file,'r') as IN,open(combined_cell_metrics_file+'.clean','w') as OUT:
        for line in IN:
            line = line.strip('\n')
            if line.startswith('Cell'):
                contents = line.split('\t')
                contents[1:6] = clean_header_metrics
                print >> OUT,"\t".join(contents)
                continue
            else:                
                contents = line.split('\t')
                cell = contents[0]
                print >> OUT,line

    with open(combined_umi_counts_file,'r') as IN,open(combined_umi_counts_file+'.clean','w') as OUT:
        i = 0
        for line in IN:
            line = line.strip('\n')
            contents = line.split('\t')
            if i == 0: ## Header
                header_anno = clean_header_umi
                header_cells = contents[6:]
                outheader = '\t'.join(header_anno) + '\t' + '\t'.join(header_cells)
                print >> OUT,outheader
                i=i+1
                continue
            else:                
                print >> OUT,line
           
def sort_by_cell(outputfile,wts):
    ''' Sort the output count files by Sample_Cells

    :param: str outputfile: path to the output file
    :param: bool wts: whether the file corresponds to a primer/gene file
    '''
    temp=outputfile+'.sorted'
    with open(outputfile,'r') as IN,open(temp,'w') as OUT:
        i=0
        for line in IN:
            contents=line.strip('\n').split('\t')
            if i == 0: #header
                contents = line.strip('\n').split('\t')
                if wts:
                    anno_header = '\t'.join(contents[0:6])
                    cells = contents[6:]
                else:
                    anno_header = '\t'.join(contents[0:7])
                    cells = contents[7:]                    
                    
                sorted_cells = '\t'.join(natsort.natsorted(cells))
                OUT.write(anno_header+'\t'+sorted_cells+'\n')
                i=1
                continue
            if wts:
                anno = '\t'.join(contents[0:6])
                sorted_vals = [x for _,x in natsort.natsorted(zip(cells,contents[6:]))]
            else:
                anno = '\t'.join(contents[0:7])
                sorted_vals = [x for _,x in natsort.natsorted(zip(cells,contents[7:]))]               
            OUT.write(anno+'\t'+'\t'.join(sorted_vals)+'\n')
    os.system('mv {temp} {outputfile}'.format(temp=temp,outputfile=outputfile))

def read_cell_file(cfile,metric_dict,is_lowinput):
    ''' Read a cell metrics file and return the parsed metrics as a dict

    :param str cfile: the cell metric file
    :param dict metric_dict: a dict of dict of metrics
    :return the dictionary of parsed metrics
    :param str: is_lowinput: Whether the protocol was for a low input application(1/0)
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
        sample = os.path.basename(metric_file).rstrip("_read_stats.txt")
        assert sample != '', "Error could not identify sample name from file path : {}".format(metric_file)
        for line in IN:
            metric,val = line.strip('\n').split(':')
            metric_dict[metric][sample] = float(val)
    return metric_dict

def check_metric_counts(sample_metrics,cell_metrics,UMI_gene_count):
    ''' Sanity check to ensure metrics tally up between sampleindex and cellindex
    level metrics
    :param dict sample_metrics: dict containing metrics aggregated over all samples
    :param dict cell_metrics: dict containing metrics aggregated over all cell indices
    :param int UMI_gene_count: total UMI count over all genes
    '''
    assert sample_metrics['reads used'] ==  cell_metrics['reads used'],"{sample} != {cell} , Read accounting failed !".format(sample=sample_metrics['reads used'],cell=cell_metrics['reads used']) 
    assert sample_metrics['total UMIs'] == cell_metrics['UMIs'],"UMI accounting failed !"
    assert sample_metrics['total UMIs'] == UMI_gene_count,"UMI accounting failed !"    

def combine_sample_metrics(files_to_merge,outfile,is_lowinput,cells_dropped,output_dir):
    ''' Combine metrics on the sample level similar to the cells
    :param list files_to_merge: the files to merge
    :param outfile: the output file to write to
    :param str: is_lowinput: Whether the protocol was for a low input application(1/0)
    :param list cells_dropped: cells which were dropped
    :param str output_dir: base output directory, use this for searching for the read statistic
                           files for the cells dropped

    :return dict containing some metrics aggregated over all samples to be used for read accounting
    :rtype dict
    '''
    sample_metrics = MyOrderedDict()
    dropped_metrics = defaultdict(lambda:defaultdict(int))
    new_metric = 'reads dropped, cell has no genes with more than 5 UMIs'    
    ## Get UMIs and reads used from the cells which were dropped
    for cell in cells_dropped:
        if len(cell.split('_')) == 2:
            sample_index,cell_index = cell.split('_')            
        else:    
            sample_index = '_'.join(cell.split('_')[0:-1])
            cell_index = cell.split('_')[-1]

        read_stats_file = glob.glob(os.path.join(output_dir,'*/{sample_index}/Cell{cell_index}_*/read_stats.txt'.format(sample_index=sample_index,cell_index=cell_index)))[0]
        cell_stats_file = glob.glob(os.path.join(output_dir,'*/{sample_index}/Cell{cell_index}_*/*_demultiplex_stats.txt'.format(sample_index=sample_index,cell_index=cell_index)))[0]
        
        cell_check = os.path.dirname(read_stats_file).split('/')[-1].split('_')[0].strip('Cell')
        ## Check to make sure we got the correct file
        assert cell_check == cell_index, "Incorrect matching of dropped CellIndex : {}".format(cell+" to "+cell_check)
        reads_dropped_less_5_UMI=0
        UMIs_dropped=0
        with open(read_stats_file,'r') as IN:
            for line in IN:                
                metric,val = line.strip('\n').split(':')
                if metric!="detected genes":
                    dropped_metrics[metric][sample_index]+=int(val)

        with open(cell_stats_file,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                if metric == 'reads total':
                    found=True
                    dropped_metrics[new_metric][sample_index]+= int(val)
                    reads_total = int(val)
                else:
                    after_qc = int(val)
        dropped_metrics['reads dropped, less than 25 bp'][sample_index]+= reads_total - after_qc
        
    ## Read metrics for each sample
    for sfile in files_to_merge:
        sample_metrics = read_sample_metrics(sfile,sample_metrics)
    ## Update sample_metrics to account for cells dropped
    total_used_reads = defaultdict(int)
    for metric in dropped_metrics:
        if metric!=new_metric:
            for sample_index in dropped_metrics[metric]:
                if metric.startswith('reads used,'):
                    total_used_reads[sample_index]+=sample_metrics[metric][sample_index]
                sample_metrics[metric][sample_index] = sample_metrics[metric][sample_index] - \
                                                       dropped_metrics[metric][sample_index]
                
    for sample_index in dropped_metrics[metric]: ## Update mean reads per UMI to account for dropping
        sample_metrics['mean reads per UMI'][sample_index] = total_used_reads[sample_index]/float(sample_metrics['total UMIs'][sample_index])
                
    ## Combine and write resultant output file
    return_metrics = defaultdict(int)
    with open(outfile,'w') as OUT:
        i=0
        for metric in sample_metrics:
            if i == 0:
                header = 'Samples\t'+'\t'.join(sample_metrics[metric].keys())
                OUT.write(header+'\n')
                i+=1
            out = metric
            for sample in sample_metrics[metric]:
                if metric.startswith('reads used,'):
                    return_metrics['reads used']+=int(sample_metrics[metric][sample])
                elif metric.startswith('total UMIs'):
                    return_metrics['total UMIs']+=int(sample_metrics[metric][sample])                   
                out = out+'\t'+float_to_string(round(sample_metrics[metric][sample],2))               
            
            OUT.write(out+'\n')
            if metric in ['reads dropped, less than 25 bp endogenous seq after primer','reads dropped, aligned to genome, multiple loci']:
                ## Add new metric for cells dropped                
                out = new_metric
                for sample in sample_metrics[metric]:
                    if sample in dropped_metrics[new_metric]:
                        out = out+'\t'+float_to_string(round(dropped_metrics[new_metric][sample],2))
                    else:
                        out = out+'\t'+float_to_string(0.0)
                OUT.write(out+'\n')                

    return return_metrics

def combine_cell_metrics(files_to_merge,outfile,is_lowinput,cells_to_restrict):
    ''' Combine cell metrics from different samples
    :param list files_to_merge: the files to merge
    :param str outfile: the outputfile to write the aggregate metrics
    :param str: is_lowinput: Whether the protocol was for a low input application(1/0)
    :param list cells_to_restrict: restrict cells to this list
    
    :return Dict containing aggregated metrics over all cells
    :rtype dict
    '''    
    cell_metrics = MyOrderedDict()
    files_to_merge  = natsort.natsorted(files_to_merge)
    for cfile in files_to_merge:
        cell_metrics = read_cell_file(cfile,cell_metrics,is_lowinput)

    return_metrics = defaultdict(int)
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
                if metric.startswith('reads used,') or metric == 'UMIs':
                    if metric.startswith('reads used,'):
                        met = 'reads used'
                    else:
                        met = 'UMIs'
                    return_metrics[met]+=int(cell_metrics[cell][metric])                    
                out = out+'\t'+cell_metrics[cell][metric]
            OUT.write(out+'\n')
            
    return return_metrics

def combine_count_files(files_to_merge,outfile,wts,cells_to_restrict=[]):
    ''' Function to combine cells from different samples into 1 file
    The directory strucuture of a typical run loooks like :
            <run_dir>
              <primary_analysis>
                --- combined_umi_count_file.txt
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
            cells which were dropped
            sum of UMI count over all cells
    :rtype tuple of (list,list,int) 
    ''' 
    i = 0
    UMI = defaultdict(lambda:defaultdict(int))
    header_cells = set()
    cells_dropped = set()
    ## Iterate over the files to merge
    for f in files_to_merge:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        sample_name = os.path.dirname(f).split('/')[-2]
        check_counts = []
        with open(f,'r') as IN:
            cell_key = sample_name+'_'+str(cell)
            for line in IN:
                if wts:
                    k1,k2,k3,k4,k5,k6,umi = line.rstrip('\n').split('\t')
                    key = (k1,k2,k3,k4,k5,k6)
                else:
                    k1,k2,k3,k4,k5,k6,k7,umi = line.rstrip('\n').split('\t')
                    key = (k1,k2,k3,k4,k5,k6,k7)
                ## Hash the umi counts by annotation and cells
                UMI[key][cell_key] = umi
                check_counts.append(int(umi))
                
            if not wts:
                if cell_key in cells_to_restrict:
                    header_cells.add(cell_key)
                else:
                    cells_dropped.add(cell_key)
            else:
                if any(e >= 5 for e in check_counts): ## Check to make sure the cell has atleast 5 UMI count for any 1 gene
                    header_cells.add(cell_key)
                else:
                    cells_dropped.add(cell_key)
    ## Create header
    if wts:
        header = "gene id\tgene\tstrand\tchrom\tloc 5' GRCh38\tloc 3' GRCh38\t{cells}\n"
    else:
        header = "gene id\tgene\tstrand\tchrom\tloc 5' GRCh38\tloc 3' GRCh38\tprimer seq\t{cells}\n"
    temp = '\t'.join(list(header_cells))
    head = header.format(cells=temp)
    ## Print output
    total_UMIs = 0
    with open(outfile,'w') as OUT:
        OUT.write(head)
        for key in UMI:
            write=True
            out = '\t'.join(key)
            for cell in header_cells:
                if cell not in UMI[key]:
                    raise Exception("Cell not hashed for Gene/Primer : {cell}-{k}".format(cell=cell,k=key))
                else:
                    out = out + '\t{}'.format(UMI[key][cell])                    
                    total_UMIs+=int(UMI[key][cell])
            if write:        
                OUT.write(out+'\n')                
    ## Sort the count file
    sort_by_cell(outfile,wts)

    return (header_cells,cells_dropped,total_UMIs)
