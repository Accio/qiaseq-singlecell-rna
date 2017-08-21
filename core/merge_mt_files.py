import glob
import os
import sys
from collections import defaultdict,OrderedDict

def merge_count_files_wts(basedir,out_file,sample_name,ncells=384):
    ''' Merge count files from different cells
    :param str basedir: the basedirectory
    :param str out_file: the path to the output file
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
                chrom,start,stop,strand,gene,gene_type,mt = line.rstrip('\n').split('\t')
                key = (chrom,start,stop,strand,gene,gene_type)
                MT[key][cell] = mt
    last_key = cell

    with open(out_file,'w') as OUT:
        print >> OUT,"chromosome\tstart\tstop\tstrand\tgene\tgene_type\t{cells}".format(cells=cell_header)
        for key in MT:
            out = '\t'.join(key)
            for cell in cells:
                if str(cell) not in MT[key]:
                    out = out + '\t0'
                else:
                    out = out + '\t'+MT[key][str(cell)]
            print >> OUT,out

def merge_metric_files_wts(basedir,temp_metric_file,metric_file,metric_file_cell,sample_name,ncells=384):
    ''' Merge the metrics from gene finding
    :param str basedir: the basedirectory
    :param str temp_metric_file: the temp metric file with the overall demultiplex stats
    :param str metric_file: the path to the metric file
    :param str metric_file_cell: the path to the metric file with info for each cell
    :param str sample_name: the sample name
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

    ## Get Alignment Stats
    files = glob.glob(os.path.join(basedir,"*/read_stats.txt"))
    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                if metric != 'num_genes_annotated':
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
        ## Calc percentage of reads for which genes were annotated
        temp=(float(metric_dict['num_reads_used'])/metric_dict['num_reads_demultiplexed_for_alignment'])*100
        M.write('perc_reads_used: {}\n'.format(float_to_string(round(temp,2))))

    ## Write metrics for each cell to file
    write_metrics_cells_wts(metric_dict_per_cell,ncells,sample_name,metric_file_cell)

def float_to_string(val):
    return ('%.2f' % val).rstrip('0').rstrip('.')

def write_metrics_cells_wts(metric_dict_per_cell,ncells,sample_name,outfile):
    ''' Write metrics for each cell

    :param dict of dict: metric_dict_per_cell: ['cell']['metric'] -> val
    :param int: ncells: the number of cells
    :param str: sample_name: the name of the sample
    :param str: outfile: the output file to write to
    :return None
    '''
    cells = range(1,ncells+1)
    with open(outfile,'w') as OUT:
        header = (
            "Cell\tTotal_reads\tTotal_pass_QC_reads\tDetected_genes\tMapped_reads\t"
            "Uniquely_Mapped_reads\tUsed_reads\tUsed_reads_ERCC\tMapped_ratio\t"
            "Used_ratio\n"
        )
        OUT.write(header)
        for cell in cells:
            cell = str(cell)
            key = 'Cell'+cell+'_'+sample_name
            if cell not in metric_dict_per_cell: ## Cell had no reads in demultiplexing
                out = key+'\t'+'\t'.join(['0']*10)+'\n'
                OUT.write(out)
            else:
                num_reads = metric_dict_per_cell[cell]['num_reads']
                after_qc_reads = metric_dict_per_cell[cell]['after_qc_reads']
                num_reads_mapped = metric_dict_per_cell[cell]['num_reads_mapped']
                num_reads_multimapped = metric_dict_per_cell[cell]['num_reads_multimapped']
                num_reads_unknown_chrom = metric_dict_per_cell[cell]['num_reads_unknown_chrom']
                num_reads_not_annotated = metric_dict_per_cell[cell]['num_reads_not_annotated']
                num_reads_used = metric_dict_per_cell[cell]['num_reads_used']
                num_reads_used_ercc = metric_dict_per_cell[cell]['num_reads_used_ercc']
                num_genes_annotated = metric_dict_per_cell[cell]['num_genes_annotated']
                out = [
                    key,str(num_reads),str(after_qc_reads),
                    str(num_genes_annotated),str(num_reads_mapped),
                    str(num_reads_mapped-num_reads_multimapped),
                    str(num_reads_used),str(num_reads_used_ercc),
                    float_to_string(round(float(num_reads_mapped)/after_qc_reads,2)),
                    float_to_string(round(float(num_reads_used)/after_qc_reads,2))
                ]
                OUT.write('\t'.join(out))
                OUT.write('\n')
def merge_count_files(basedir,out_file,sample_name,ncells=96):
    '''
    :param str basedir: the basedirectory
    :param str out_file: the path to the output file
    :param int ncells: the number of cell indices
    '''
    i = 0
    MT = defaultdict(lambda:defaultdict(list))
    cells = range(1,ncells+1)
    cell_header = '\t'.join(sample_name+'_Cell'+str(cell) for cell in cells)
    files = glob.glob(os.path.join(basedir,"*/mt_count.txt"))

    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                chrom,start,stop,seq,strand,gene,mt = line.rstrip('\n').split('\t')
                MT[seq][cell] = [chrom,start,stop,seq,gene,mt]

    last_key = cell ## Store the last cell from the previous loop to get primer info from the dictionary
    with open(out_file,'w') as OUT:
        print >> OUT,"chromosome\tstart\tstop\tsequence\tgene\t%s"%cell_header
        for key in MT:
            out='\t'.join(MT[key][last_key][0:-1])
            for cell in cells:
                if str(cell) not in MT[key]:
                    out = out + '\t0'
                else:
                    #print cell
                    out=out+'\t'+MT[key][str(cell)][-1]
            print >> OUT,out

def merge_metric_files(basedir,metric_file,ncells=96):
    ''' Merge the primer finding metrics
    :param str basedir: the basedirectory
    :param str metic_file: the path to the metric file
    :param int ncells: the number of cell indices
    '''
    i = 0
    cells = range(1,ncells+1)
    files = glob.glob(os.path.join(basedir,"*/read_stats.txt"))
    metric_dict = defaultdict(int)

    for f in files:
        cell = os.path.dirname(f).split('/')[-1].split('_')[0].strip('Cell')
        with open(f,'r') as IN:
            for line in IN:
                metric,val = line.strip('\n').split(':')
                metric_dict[metric]+=int(val)

    temp_dict = {}
    with open(metric_file,'r') as M:
        for line in M:
            metric,val = line.strip('\n').split(':')
            temp_dict[metric] = float(val)

    with open(metric_file,'a') as M:
        for metric,val in metric_dict.items():
            print >> M,metric+': '+str(val)

        ## Calc percentage of reads for which primers were found
        temp=(float(metric_dict['num_reads_primer_found'])/temp_dict['num_reads'])*100
        print >> M,'perc_reads_found_primer: '+str(temp)

if __name__ == '__main__':
    merge_count_files(sys.argv[1],sys.argv[2])
