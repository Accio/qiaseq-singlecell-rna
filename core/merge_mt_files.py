import glob
import os
import sys
from collections import defaultdict,OrderedDict

        
def float_to_string(val):
    ''' Helper function to convert a float into a printable string

    :param float val: the value to typecast
    '''
    return ('%.2f' % val).rstrip('0').rstrip('.')

def merge_count_files(basedir,out_file,sample_name,wts,ncells=384):
    ''' Merge count files from different cells

    :param str basedir: the basedirectory
    :param str out_file: the path to the output file
    :param str sample_name: the name of the sample
    :param bool wts: Whether it is whole transcriptome or not
    :param int ncells: the number of cell indices
    '''
    i = 0
    MT = defaultdict(lambda:defaultdict(int))
    cells = range(1,ncells+1)
    cell_header = '\t'.join(sample_name+'_Cell'+str(cell) for cell in cells)
    files = glob.glob(os.path.join(basedir,"*/mt_count.txt"))


    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                k1,k2,k3,k4,k5,k6,mt = line.rstrip('\n').split('\t')
                key = (k1,k2,k3,k4,k5,k6)
                MT[key][cell] = mt
    last_key = cell

    if wts:
        header = "chromosome\tstart\tstop\tstrand\tgene\tgene_type\t{cells}\n"
    else:
        header = "chromosome\tstart\tstop\tstrand\tgene\tprimer_sequence\t{cells}\n"
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

def merge_metric_files(basedir,temp_metric_file,metric_file,metric_file_cell,sample_name,wts=True,ncells=384):
    ''' Merge the metrics from primer/gene finding

    :param str basedir: the basedirectory
    :param str temp_metric_file: the temp metric file with the overall demultiplex stats
    :param str metric_file: the path to the metric file
    :param str metric_file_cell: the path to the metric file with info for each cell
    :param str sample_name: the sample name
    :param bool: wts: Whether whole transcriptome or not
    :param int ncells: the number of cell indices
    '''
    i = 0
    cells = range(1,ncells+1)
    metric_dict = OrderedDict()
    metric_dict_per_cell = defaultdict(lambda:defaultdict(int))

    ## Read temp metrics file
    with open(temp_metric_file,'r') as M:
        for line in M:
            metric,val = line.strip('\n').split(':')
            metric_dict[metric] = float(val)

    choice = 'num_genes_annotated' if wts else 'num_primers_found'
    
    ## Get Alignment Stats
    files = glob.glob(os.path.join(basedir,"*/read_stats.txt"))
    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                if metric != choice:
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


    with open(metric_file,'w') as M:
        for metric,val in metric_dict.items():
            out = metric+': {}\n'.format(float_to_string(round(val,2)))
            M.write(out)

        choice = 'num_reads_used' if wts else 'num_reads_primer_found'
        ## Calc percentage of reads for which genes/primers were annotated
        temp=(float(metric_dict[choice])/metric_dict['num_reads_demultiplexed_for_alignment'])*100
        M.write('perc_reads_used: {}\n'.format(float_to_string(round(temp,2))))

    ## Write metrics for each cell to file
    write_metrics_cells(metric_dict_per_cell,ncells,sample_name,metric_file_cell,wts)

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
            "Cell\tTotal_reads\tTotal_pass_QC_reads\tDetected_genes\t"
            "Mapped_reads\tUniquely_Mapped_reads\tUsed_reads\tUsed_reads_ERCC"
            "\tMapped_ratio\tUsed_ratio\n"
        )
        header_len = len(header.split('\t'))
    else:            
        header = (
            "Cell\tTotal_reads\tTotal_pass_QC_reads\tDetected_primers\t"
            "Primers_offtarget\tPrimers_Mismatch\tPrimers_offloci\t"
            "Primers_unmapped\tUsed_ratio\n"
        )
        header_len = len('\t'.join(header)) 
    with open(outfile,'w') as OUT:
        OUT.write(header)
        for cell in cells:
            cell = str(cell)
            key = 'Cell'+cell+'_'+sample_name
            if cell not in cell_metrics: ## Cell had no reads in demultiplexing
                out = key+'\t'+'\t'.join(['0']*(header_len-1))+'\n'
                OUT.write(out)
            else:
                if wts:
                    out = [key , str(cell_metrics[cell]['num_reads']),
                        str(cell_metrics[cell]['after_qc_reads']),
                        str(cell_metrics[cell]['num_genes_annotated']),
                        str(cell_metrics[cell]['num_reads_mapped']),
                        str(cell_metrics[cell]['num_reads_mapped'] - \
                            cell_metrics[cell]['num_reads_multimapped']),
                        str(cell_metrics[cell]['num_reads_used']),
                        str(cell_metrics[cell]['num_reads_used_ercc']),
                        float_to_string(
                            round(float(
                                cell_metrics[cell]['num_reads_mapped']) / \
                                  cell_metrics[cell]['after_qc_reads'],2)),
                        float_to_string(
                            round(float(
                                cell_metrics[cell]['num_reads_used']) / \
                                  cell_metrics[cell]['after_qc_reads'],2))
                    ]
                else:
                    out = [ key, str(cell_metrics[cell]['num_reads']),
                        str(cell_metrics[cell]['after_qc_reads']),
                        str(cell_metrics['num_primers_found']),
                        str(cell_metrics['num_reads_primer_offtarget']),
                        str(cell_metrics['num_reads_primer_mismatch']),
                        str(cell_metrics['num_reads_primer_off_loci']),
                        str(cell_metrics['num_reads_unmapped']),
                        str(cell_metrics['num_reads_endogenous_seq_not_matched']),
                        str(cell_metrics['num_reads_primer_found']),
                        float_to_string(
                            round(float(
                                cell_metrics[cell]['num_reads_primer_found']) / \
                                  cell_metrics[cell]['after_qc_reads'],2))
                    ]

                    
                OUT.write('\t'.join(out))
                OUT.write('\n')
                assert header_len == len(out), "Error in Column Lengths!!"
                
if __name__ == '__main__':
    merge_count_files(sys.argv[1],sys.argv[2])