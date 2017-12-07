import glob
import os
import sys
from collections import defaultdict,OrderedDict

        
def float_to_string(val):
    ''' Helper function to convert a float into a printable string

    :param float val: the value to typecast
    '''
    return ('%.2f' % val).rstrip('0').rstrip('.')

def merge_count_files(basedir,out_file,sample_name,wts,ncells,files_to_merge):
    ''' Merge count files from different cells

    :param str basedir: the basedirectory
    :param str out_file: the path to the output file
    :param str sample_name: the name of the sample
    :param bool wts: Whether it is whole transcriptome or not
    :param int ncells: the number of cell indices
    :param list files_to_merge: which files to merge
    '''
    i = 0
    MT = defaultdict(lambda:defaultdict(int))
    cells = range(1,ncells+1)
    cell_header = '\t'.join(sample_name+'_Cell'+str(cell) for cell in cells)
    for f in files_to_merge:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                k1,k2,k3,k4,k5,k6,mt = line.rstrip('\n').split('\t')
                key = (k1,k2,k3,k4,k5,k6)
                MT[key][cell] = mt
    last_key = cell

    if wts:
        header = "gene id\tgene\tstrand\tchrom\tloc 5' GRCH38\tloc 3' GRCH38\t{cells}\n"
    else:
        header = "primer sequence\tgene\tstrand\tchrom\tloc 5' GRCH38\tloc 3' GRCH38\t{cells}\n"
    with open(out_file,'w') as OUT:
        OUT.write(header.format(cells = cell_header))
        for key in MT:
            out = '\t'.join(key)
            for cell in cells:
                if str(cell) not in MT[key]:
                    out = out + '\t0'
                else:
                    out = out + '\t'+MT[key][str(cell)]
            OUT.write(out+'\n')

def merge_metric_files(basedir,temp_metric_file,metric_file,metric_file_cell,sample_name,wts,ncells,files_to_merge):
    ''' Merge the metrics from primer/gene finding

    :param str basedir: the basedirectory
    :param str temp_metric_file: the temp metric file with the overall demultiplex stats
    :param str metric_file: the path to the metric file
    :param str metric_file_cell: the path to the metric file with info for each cell
    :param str sample_name: the sample name
    :param bool: wts: Whether whole transcriptome or not
    :param int ncells: the number of cell indices
    :param list files_to_merge: which files to merge
    '''
    i = 0
    cells = range(1,ncells+1)
    metric_dict = OrderedDict()
    metric_dict_per_cell = defaultdict(lambda:defaultdict(int))
    do_not_add_metrics = ['detected genes']
    ## Read temp metrics file to get total reads , reads dropped during demultiplexing
    with open(temp_metric_file,'r') as M:
        for line in M:
            metric,val = line.strip('\n').split(':')
            metric_dict[metric] = float(val)
    ## Get read stats for each cell as well as aggregating for sample level
    for f in files_to_merge:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                if metric not in do_not_add_metrics:
                    if metric not in metric_dict:
                        metric_dict[metric]=int(val)
                    else:
                        metric_dict[metric]+=int(val)
                metric_dict_per_cell[cell][metric] = int(val)
    ## Get Per Cell Demultiplex Stats
    files = glob.glob(os.path.join(basedir,"*/*_demultiplex_stats.txt"))
    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                metric_dict_per_cell[cell][metric] = int(val)

    ## Write metrics for the sample
    write_metrics_sample(metric_dict,metric_file,wts)
    ## Write metrics for each cell to a file
    write_metrics_cells(metric_dict_per_cell,ncells,sample_name,metric_file_cell,wts)
    
def write_metrics_sample(sample_metrics,outfile,wts):
    ''' Write metrics on sample level

    :param dict sample_metrics: <metric> -> <val>
    :param str outfile: the outputfile to write the metrics
    :param bool wts: whether this is whole transcriptome seq
    '''
    ## Check to make sure the metrics add up
    reads_total = int(sample_metrics['reads total'])    
    reads_dropped_demultiplexing = (
        int(sample_metrics['reads dropped, cell id not extracted']) + \
        int(sample_metrics['reads dropped, cell id not matching oligo']) + \
        int(sample_metrics['reads dropped, less than 25 b.p'])       
    )
    if wts:
        reads_dropped_counting = (
            int(sample_metrics['reads dropped, not mapped to genome']) + \
            int(sample_metrics['reads dropped, not annotated']) + \
            int(sample_metrics['reads dropped, aligned to genome, multiple loci']) + \
            int(sample_metrics['reads dropped, aligned to ERCC, multiple loci'])            
        )
        reads_used = (
            int(sample_metrics['reads used, aligned to genome, unique loci']) + \
            int(sample_metrics['reads used, aligned to ERCC, unique loci'])
        )        
    else:
        reads_dropped_counting = (
            int(sample_metrics['reads dropped, not mapped to genome']) + \
            int(sample_metrics['reads dropped, off target']) + \
            int(sample_metrics['reads dropped, primer not identified at read start']) + \
            int(sample_metrics['reads dropped, less than 25 b.p endogenous seq after primer'])
        )
        reads_used = (
            int(sample_metrics['reads used, aligned to genome, multiple loci']) + \
            int(sample_metrics['reads used, aligned to genome, unique loci']) + \
            int(sample_metrics['reads used, aligned to ERCC, multiple loci']) + \
            int(sample_metrics['reads used, aligned to ERCC, unique loci'])
        )        

    ## Write Overall Sample level aggregated metrics
    with open(outfile,'w') as M:
        for metric,val in sample_metrics.items():
            out = metric+': {}\n'.format(float_to_string(round(val,2)))
            M.write(out)
        ## Calc reads per UMI
        if wts:
            total_reads = float(
                int(sample_metrics['reads used, aligned to genome, unique loci']) + \
                int(sample_metrics['reads used, aligned to ERCC, unique loci'])
                )
        else:
            total_reads = float(
                int(sample_metrics['reads used, aligned to genome, unique loci']) + \
                int(sample_metrics['reads used, aligned to genome, multiple loci']) + \
                int(sample_metrics['reads used, aligned to ERCC, unique loci']) + \
                int(sample_metrics['reads used, aligned to ERCC, multiple loci'])
                )            
                
        temp=(total_reads/sample_metrics['total UMIs'])
        M.write('mean reads per UMI: {}\n'.format(float_to_string(round(temp,2))))
        
    assert (reads_total == reads_used + reads_dropped_demultiplexing + reads_dropped_counting),"Read accounting failed !"
        
def write_metrics_cells(cell_metrics,ncells,sample_name,outfile,wts):
    ''' Write metrics for each cell

    :param dict of dict: metric_dict_per_cell: ['cell']['metric'] -> val
    :param int: ncells: the number of cells
    :param str: sample_name: the name of the sample
    :param str: outfile: the output file to write to
    :param bool: wts: Whether whole transcriptome or not
    :return None
    '''
    cells = range(1,ncells+1)
    if wts:
        header = (
            "cell\treads total\treads used, aligned to genome\treads used, aligned to ERCC\tUMIs\tdetected genes\n"
        )
        header_len = len(header.split('\t'))
    else:            
        header = (
            "cell\treads total\treads used, aligned to genome\treads used, aligned to ERCC\tUMIs\tdetected genes\n"
        )
        header_len = len(header.split('\t'))
    with open(outfile,'w') as OUT:
        OUT.write(header)
        for cell in cells:
            cell = str(cell)
            key = sample_name+'_'+cell
            if cell not in cell_metrics or cell_metrics[cell]['after_qc_reads'] == 0: ## Cell had no reads in demultiplexing
                out = key+'\t'+'\t'.join(['0']*(header_len-1))+'\n'
                OUT.write(out)
            else:
                if wts:
                    reads_used_genome = str(
                        int(cell_metrics[cell]['reads used, aligned to genome, unique loci'])
                    )
                    reads_used_ercc = str(
                        int(cell_metrics[cell]['reads used, aligned to ERCC, unique loci'])
                    )
                else:
                    reads_used_genome = str(
                        int(cell_metrics[cell]['reads used, aligned to genome, unique loci']) + \
                        int(cell_metrics[cell]['reads used, aligned to genome, multiple loci'])
                    )                   
                    reads_used_ercc = str(
                        int(cell_metrics[cell]['reads used, aligned to ERCC, unique loci']) + \
                        int(cell_metrics[cell]['reads used, aligned to ERCC, multiple loci'])
                    )                    
                out = [
                    key , str(cell_metrics[cell]['reads total']),                        
                    reads_used_genome,
                    reads_used_ercc,
                    str(cell_metrics[cell]['total UMIs']),
                    str(cell_metrics[cell]['detected genes'])
                ]
                assert header_len == len(out), "Error in Column Lengths!!"
                OUT.write('\t'.join(out))
                OUT.write('\n')               
