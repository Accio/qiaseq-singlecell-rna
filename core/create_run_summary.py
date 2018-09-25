import os
import numpy as np
from xlsxwriter.workbook import Workbook
from collections import defaultdict
import ConfigParser
import glob
import gzip
import subprocess
from merge_mt_files import float_to_string

def read_sample_metrics(sample_metrics_file):
    '''
    :param str sample_metrics_file: The sample metrics file
    '''
    aggregated_metrics = dict()
    with open(sample_metrics_file,'r') as IN:
        for line in IN:
            if line.startswith('Samples') or line.startswith('mean reads per UMI'):
                continue
            contents = line.strip('\n').split('\t')
            metric = contents[0]
            values = [int(e) for e  in contents[1:]]
            total_number = sum(values)
            aggregated_metrics[metric] = total_number
            
    return aggregated_metrics

def read_cell_metrics(cell_metrics_file):
    '''
    :param str cell_metrics
    '''
    cell_metrics = defaultdict(dict)
    with open(cell_metrics_file,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            if line.startswith('Cell'):
                metrics = contents[1:]
                continue
            contents = line.strip('\n').split('\t')
            cell = contents[0]
            vals = contents[1:]
            for i in range(len(metrics)):
                cell_metrics[cell][metrics[i]] = int(vals[i])
                
    return cell_metrics

def calc_stats_gene_count(combined_gene_count_file):
    '''
    '''
    cell_metrics = defaultdict(lambda:defaultdict(int))
    temp = defaultdict(int)
    
    with open(combined_gene_count_file,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')        
            if line.startswith('gene id'):        
                cells = contents[6:]
                continue

            if contents[1].startswith('ERCC-'):
                met1 = 'num_ercc'
                met2 = 'umis_ercc'
            else:
                met1 = 'num_genes'
                met2 = 'umis_genes'

            if any(int(e) > 0 for e in contents[6:]):
                temp[met1] += 1
                temp[met2] += sum(int(e) for e in contents[6:])
                vals = contents[6:]
                for i in range(len(cells)):
                    if int(vals[i]) != 0:
                        cell_metrics[cells[i]][met1]+= 1                
                    cell_metrics[cells[i]][met2]+= int(vals[i])
                    cell_metrics[cells[i]]['umis']+= int(vals[i])

    return (cell_metrics, temp['num_genes'], temp['num_ercc'], temp['umis_genes'], temp['umis_ercc'])

def calc_median_cell_metrics(cell_metrics,metric,cells_to_drop=[],drop_outlier_cells = False):
    ''' Calculate median across all cells for a given metric
    :param dict cell_metrics: dictionary of metrics for each cell
    :param str metric: the metric to calculate the median for
    :param list: cells_to_drop: account for cell droppings, i.e. which cells to 
                                ignore when calculating the median
    :param drop_outlier_cells: whether to prune outlier cells when calculating the median
    '''
    
    ''' 
    If specified, drop some outlier cells when calculating the median
    Follow similar logic as that used in the secondary analysis, i.e.
    1.) Cells with endogenous gene UMIs below 5th percentile OR
    2.) Cells with % of endogenous UMIs below 5th percentile OR
    3.) Cells with detected genes below 5th percentile    
    '''
    if drop_outlier_cells:
        temp1 = []
        temp2 = []
        temp3 = []
        temp_cells = []
        for cell in cell_metrics:
            temp1.append(cell_metrics[cell]['umis_genes'])
            temp2.append(cell_metrics[cell]['num_genes'])
            temp3.append(float(cell_metrics[cell]['num_genes'])/(cell_metrics[cell]['num_genes'] + cell_metrics[cell]['num_ercc']))
            temp_cells.append(cell)

        num_cells = len(temp_cells)
        # compute index corresponding to the 5th percentile metrics above ; add < 5th percentile values to list of cells to drop
        for l in [temp1,temp2,temp3]:
            sorted_l = sorted(l)
            sorted_cells = [x for _,x in sorted(zip(sorted_l,temp_cells))]
            idx = int(round(0.05 * (num_cells - 1)))
            cells_to_drop.extend(sorted_cells[0:idx])
        
    temp = []
    cells_to_drop = set(cells_to_drop)
    if drop_outlier_cells:
        print "Cells dropped when computing median :\n"
        print ",".join(list(cells_to_drop))+"\n"
    
    for cell in cell_metrics:
        if cell not in cells_to_drop:
            temp.append(cell_metrics[cell][metric])

    return int(np.median(temp))
    
def is_file_empty(infile):
    ''' Helper function to check if a gzip/regular file is empty
    :param str infile: the input file
    :returns Whether the file is empty or not
    :rtype bool
    '''
    if not os.path.exists(infile):
        return True
    elif infile.endswith('.gz'):
        IN = gzip.open(infile,'r')
    else:
        IN = open(infile,'r')
    i=0
    for line in IN:
        if i == 5:
            break
        i+=1
    IN.close()
    return i == 0

def get_cells_demultiplexed(run_id):
    '''
    '''
    counter=0
    search_path = os.path.join("/home/qiauser/{run_id}/primary_analysis/*/*/*.fastq".format(run_id=run_id))
    for fastq in glob.glob(search_path):
        if not is_file_empty(fastq):
            counter+=1
    return counter


def read_clustering_cells_dropped(cells_dropped_file):
    '''
    :param str cells_dropped_file: the file path for the cells dropped from the clustering script
    :returns a list of cells dropped
    :rtype list
    '''
    ret = []
    with open(cells_dropped_file,'r') as IN:
        for line in IN:
            contents=line.strip('\n').split(',')
            if line.startswith('Cells'):
                continue
            ret.append(contents[0])
    return ret

def get_sample_names(samples_cfg):
    ''' Return sample names

    :param str samples_cfg: the samples cfg file
    :returns a comma seperated string of sample names
    :rtype str
    '''
    parser = ConfigParser.ConfigParser()
    parser.read(samples_cfg)
    ret = []
    for section in parser.sections():
        ret.append(section)
    return ','.join(ret)


def identify_DE_method_used(runid):
    ''' Identify DE method used based by grepping the log
    we use either SCDE or edgeR
    
    :param str runid: the runid corresponding to the job    
    :returns SCDE or edgeR
    :rtype str
    '''
    logfile = "/srv/qgen/code/qiaseq-singlecell-rna/{runid}.log".format(runid=runid)    
    grep_cmd = """ grep "scde based modelling failed, using edgeR" {log} """.format(log=logfile)
    try:
        out=subprocess.check_output(grep_cmd,shell=True).strip('\n')
        return 'edgeR'
    except subprocess.CalledProcessError:
        return 'SCDE'    
    
def write_run_summary(output_excel,has_clustering_run,run_id,seqtype,species,
                      samples_cfg,sample_metrics_file,cell_metrics_file,
                      cells_dropped_file,metrics_from_countfile,
                      normalization_method,hvg_method):
    '''
    :param str output_excel: the output excel file path
    :param bool has_clustering_run: whether clustering was performed for this run
    :param str run_id: the gce vm job_id
    :param str seqtype: the library/protocol type, i.e. targeted, poly-A transcriptome, drop-seq 
    :param str species: the species name, i.e. Human, Mouse, Rat
    :param str samples_cfg: the config file with sample level info like sample name
    :param str sample_metrics_file: the sample metrics file
    :param str cell_metrics_file: the cell metrics file
    :param str cells_dropped_file: file with cells dropped in the clustering script
    :param tuple metrics_from_countfile: metric stats obtained from the gene count matrix
    :param str normalization_method: the normalization scheme used in the clustering script
    :param str hvg_method: the highly variable gene selection method used    
    '''
    cell_metrics_countfile,num_genes,num_ercc,num_umis_genes,num_umis_ercc = metrics_from_countfile
    

    samples = get_sample_names(samples_cfg)
    
    if seqtype.upper() == 'WTS':
        seqtype = 'poly-A-transcriptome'
    else:
        seqtype = 'targeted'
    

    if species.upper() == 'HUMAN':
        gencode = "Gencode Release 28"
        genome = "GRCh38"
    elif species.upper() == "MOUSE":
        gencode = "Gencode Release M15"
        genome = "GRCm38"
    elif species.upper() == "RAT":
        gencode = "Ensembl Release 91"
        genome = "Ensembl Release 91"
    elif species.upper() == "MACAQUE":
        gencode = "Ensembl Release 93"
        genome = "Ensembl Release 93"
    elif species.upper() == "YEAST":
        gencode = "Ensembl Release 93"
        genome = "Ensembl Release 93"        
    else:
        raise Exception("Unsupported species ! {}".format(species))
    
    aggregated_metrics = read_sample_metrics(sample_metrics_file)
    cell_metrics = read_cell_metrics(cell_metrics_file)
    cells_demultiplexed = get_cells_demultiplexed(run_id)
    
    reads_total = aggregated_metrics['reads total']
    reads_used = 0
    for met in aggregated_metrics:
        if met.startswith('reads used'):
            reads_used += aggregated_metrics[met]

    umis = aggregated_metrics['total UMIs']
    umis_endo = num_umis_genes
    umis_ercc_controls = num_umis_ercc
    
    cells_used = len(cell_metrics.keys())
    cells_to_drop = []

    ## Account for cell droppings in the clustering step, if required
    if has_clustering_run:
        cells_to_drop = read_clustering_cells_dropped(cells_dropped_file)
        cells_used = cells_used - len(cells_to_drop)
        for cell in cells_to_drop:
            for met in cell_metrics[cell]:
                if met.startswith('reads used'):
                    reads_used = reads_used - cell_metrics[cell][met]
                    
            for met in cell_metrics_countfile[cell]:
                if met == 'umis':
                    umis = umis - cell_metrics_countfile[cell]['umis']
                elif met == 'umis_ercc':
                    umis_ercc_controls = umis_ercc_controls - cell_metrics_countfile[cell]['umis_ercc']
                elif met == 'umis_genes':
                    umis_endo = umis_endo - cell_metrics_countfile[cell]['umis_genes']
                elif met == 'num_genes':
                    num_genes = num_genes
                elif met == 'num_ercc':
                    num_ercc = num_ercc
                    
    ## Compute median values
    median_reads_per_cell = calc_median_cell_metrics(cell_metrics,'reads total',cells_to_drop)
    median_umis_per_cell = calc_median_cell_metrics(cell_metrics,'UMIs',cells_to_drop)
    median_genes_per_cell = calc_median_cell_metrics(cell_metrics_countfile,'num_genes',cells_to_drop)

    ## Open the excel file for writing and format the columns
    workbook = Workbook(output_excel)
    bold = workbook.add_format({'bold': True})
    num_fmt = workbook.add_format({'num_format': '#,###'})    
    worksheet = workbook.add_worksheet()
    ## Hard coding column width, I highly doubt this will break.
    ## PCA and K-means clustering or Gencode Release 23 are the longest words in the second column
    ## Total ERCC spike-ins detected, all cells is the longest word in the first column
    ## None of the read/umi metrics should ever exceed the length the of these two words
    ## For e.g. 10,000,000,000 (length 14 including commas is still less than 18 , length of 'Gencode Release 23')
    worksheet.set_column('A:A',len('Total ERCC spike-ins detected, all cells')+2)
    if has_clustering_run:
        worksheet.set_column('B:B',len('PCA and K-means clustering')+2)
    else:        
        worksheet.set_column('B:B',len('Gencode Release 23')+2)

    ## Write the summary data    
    worksheet.write("A1","Run level summary",bold)    
    worksheet.write_row(1,0,["Job identifier",run_id])
    worksheet.write_row(2,0,["Library protocol",seqtype])    
    worksheet.write_row(3,0,["Reference genome",genome])
    worksheet.write_row(4,0,["Transcriptome models",gencode])
    worksheet.write_row(5,0,["Read mapping","STAR version 2.5.3a"])
    
    worksheet.write_row(6,0,[]) ## Blank    
    worksheet.write(7,0,"UMI and read counts",bold)
    
    worksheet.write(8,0,"Read fragments, total")
    worksheet.write_number(8,1,reads_total,num_fmt)
    
    worksheet.write(9,0,"Read fragments, used")
    worksheet.write_number(9,1,reads_used,num_fmt)
    
    worksheet.write_row(10,0,["Read fragments per UMI",float(float_to_string(float(reads_used)/umis))])
    
    worksheet.write(11,0,"UMIs")
    worksheet.write_number(11,1,umis,num_fmt)
    
    worksheet.write(12,0,"UMIs, endogenous genes")
    worksheet.write_number(12,1,umis_endo,num_fmt)
    
    worksheet.write(13,0,"UMIs, ERCC controls")
    worksheet.write(13,1,umis_ercc_controls,num_fmt)
    
    worksheet.write_row(14,0,[]) ## Blank    
    worksheet.write(15,0,"Cell and gene level summary",bold)
    
    worksheet.write(16,0,"Cells demultiplexed")
    worksheet.write_number(16,1,cells_demultiplexed,num_fmt)
    
    worksheet.write(17,0,"Cells used")
    worksheet.write_number(17,1,cells_used,num_fmt)
    
    worksheet.write(18,0,"Median read fragments per cell")
    worksheet.write_number(18,1,median_reads_per_cell,num_fmt)
    
    worksheet.write(19,0,"Median UMIs per cell")
    worksheet.write_number(19,1,median_umis_per_cell,num_fmt)
    
    worksheet.write(20,0,"Median genes per cell")
    worksheet.write_number(20,1,median_genes_per_cell,num_fmt)
    
    worksheet.write(21,0,"Total genes detected, all cells")
    worksheet.write_number(21,1,num_genes,num_fmt)
    
    worksheet.write(22,0,"Total ERCC spike-ins detected, all cells")
    worksheet.write_number(22,1,num_ercc,num_fmt)
    
    worksheet.write_row(23,0,[]) ## Blank
    if has_clustering_run:
        worksheet.write(24,0,"Expression analysis",bold)
        worksheet.write_row(25,0,["UMI count normalization",normalization_method])
        worksheet.write_row(26,0,["Highly variable gene selection",hvg_method])
        worksheet.write_row(27,0,["Cell clustering","PCA and K-means clustering"])
        worksheet.write_row(28,0,["Differential expression analysis",identify_DE_method_used(run_id)])
    workbook.close()
