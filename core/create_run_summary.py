import os
import numpy as np
from xlsxwriter.workbook import Workbook
from collections import defaultdict
import ConfigParser
import glob
import gzip

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
            cell_id = contents[0]
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

            if contents[1].startswith('ERCC'):
                met1 = 'num_ercc'
                met2 = 'umis_ercc'
            else:
                met1 = 'num_genes'
                met2 = 'umis_genes'

            if any(int(e) > 0 for e in contents[6:]):
                temp[met1]+=1
                temp[met2] = sum(int(e) for e in contents[6:])
                vals = contents[6:]
                for i in range(len(cells)):
                    if int(val[i]) != 0:
                        cell_metrics[cells[i]][met1]+= 1                
                    cell_metrics[cells[i]][met2]+= int(val[i])
                    cell_metrics[cells[i]]['UMIs']+= int(val[i])
                
    return (cell_metrics, temp['num_genes'], temp['num_ercc'], temp['umis_genes'], temp['umis_ercc'])

def calc_median_cell_metrics(cell_metrics,metric,cells_to_drop=[]):
    ''' Calculate median across all cells for a given metric
    :param dict cell_metrics: dictionary of metrics for each cell
    :param str metric: the metric to calculate the median for
    :param list: cells_to_drop: account for cell droppings, i.e. which cells to 
                                ignore when calculating the median
    '''
    for cell in cell_metrics:
        if cell not in cells_to_drop:
            temp.append(cell_metrics[cell][metric])

    return int(np.median(temp))
    
def is_gzip_empty(gzipfile):
    ''' Helper function to check if a gzip file is empty
    :param str gzipfile: the .gz file
    :returns Whether the file is empty or not
    :rtype bool
    '''
    with gzip.open(gzipfile,'r') as IN:
        i=0
        for line in IN:
            if i == 5:
                break
            i+=1
    return i == 0

def get_cells_demultiplexed(run_id):
    '''
    '''
    counter=0
    search_path = os.path.join("/home/qiauser/{run_id}/primary_analysis/*/*/*.fastq.gz")
    for fastq in glob.glob(search_path):
        if not is_gzip_empty(fastq):
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
    parser.read(self.samples_cfg)
    ret = []
    for section in parser.sections():
        ret.append(section)
    return ','.join(ret)

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
    
    workbook = Workbook(output_excel)
    bold = workbook.add_format({'bold': True})    
    worksheet = workbook.add_worksheet()

    samples = get_sample_names(samples_cfg)
    
    worksheet.write("A1","Run level summary",bold)
    worksheet.write_row(1,0,["Job identifier",run_id])
    worksheet.write_row(2,0,["Read set names",samples])
    worksheet.write_row(3,0,["Library protocol",seqtype])

    if species.upper() == 'HUMAN':
        gencode = "Gencode Release 23"
        genome = "GRCh38"
    elif species.upper() == "MOUSE":
        gencode = "Gencode Release M15"
        genome = "GRCm38"
    else:
        raise Exception("Unsupported species ! {}".format(species))
    
    worksheet.write_row(4,0,["Reference genome",genome])
    worksheet.write_row(5,0,["Transcriptome models",gencode])
    worksheet.write_row(6,0,["Read mapping","STAR version 2.5.3a"])
    
    worksheet.write_row(7,0,[]) ## Blank

    aggregated_metrics = read_sample_metrics(sample_metrics_file)
    cell_metrics = read_cell_metrics(cell_metrics_file)
    cells_demultiplexed = get_cells_demultiplexed(run_id)
    
    reads_total = aggregated_metrics['reads total']/2
    for met in aggregated_metrics:
        if met.startswith('reads used'):
            reads_used += aggregated_metrics[met]

    umis = aggregated_metrics['total UMIs']
    umis_endo = num_umis_genes
    umis_ercc_controls = num_umis_ercc
    cells_used = len(cell_metrics.keys())
    cells_to_drop = []
    ## Account for cell droppings in clustering stage if required
    if has_clustering_run:
        cells_to_drop = read_clustering_cells_dropped(cells_dropped_file)
        cells_used = cells_used - len(cells_to_drop)
        for cell in cells_to_drop:
            for met in cell_metrics[cell]:
                if met.startswith('reads total'):
                    reads_total = reads_total - cell_metrics[cell][met]
                if met.startswith('reads used'):
                    reads_used = reads_used - cell_metrics[cell][met]

            for met in cell_metrics_countfile[cell]:
                if met == 'umis':
                    umis = umis_ercc_controls - cell_metrics_countfile['umis']
                if met == 'umis_ercc':
                    umis_ercc_controls = umis_ercc_controls - cell_metrics_countfile['umis_ercc']
                if met == 'umis_genes':
                    umis_endo = umis_endo - cell_metrics_countfile['umis_genes']
                if met == 'num_genes':
                    num_genes = num_genes - cell_metrics_countfile['num_genes']
                if met == 'num_ercc':
                    num_ercc = num_genes - cell_metrics_countfile['num_ercc']

    ## Compute median values
    median_reads_per_cell = calc_median_cell_metrics(cell_metrics,'reads total',cells_to_drop)
    median_umis_per_cell = calc_median_cell_metrics(cell_metrics,'UMIs',cells_to_drop)
    median_genes_per_cell = calc_median_cell_metrics(cell_metrics_countfile,'num_genes',cells_to_drop)
    
    
    worksheet.write("A9","UMI and read counts",bold)
    worksheet.write("A10","Read fragments, total")
    worksheet.write("B10",reads_total)
    worksheet.write("A11","Read fragments, used")
    worksheet.write("B11",reads_used)
    worksheet.write("A12","Read fragments per UMI")
    worksheet.write("B12",float(float_to_string(float(reads_used)/umis)))
    worksheet.write("A13","UMIs")
    worksheet.write("B13",umis)
    worksheet.write("A14","UMIs, endogenous genes")
    worksheet.write("B14",umis_endo)
    worksheet.write("A15","UMIs, ERCC controls")
    worksheet.write("B15",umis_ercc_controls)    
    
    worksheet.write_row(15,0,[]) ## Blank

    worksheet.write("A17","Cell and gene level summary",bold)
    worksheet.write("A18","Cells demultiplexed")
    worksheet.write("B18",cells_demultiplexed)
    worksheet.write("A19","Cells used")
    worksheet.write("B19",cells_used)
    worksheet.write("A20","Median read fragments per cell")
    worksheet.write("B20",median_reads_per_cell)
    worksheet.write("A21","Median UMIs per cell")
    worksheet.write("B21",median_umis_per_cell)
    worksheet.write("A22","Median genes per cell")
    worksheet.write("B22",median_genes_per_cell)
    worksheet.write("A23","Total genes detected, all cells")
    worksheet.write("B23",num_genes)
    worksheet.write("A23","Total ERCC spike-ins detected, all cells")
    worksheet.write("B23",num_ercc)    

    worksheet.write_row(23,0,[])
    
    if not is_low_input:
        worksheet.write("A26","Expression analysis",bold)
        worksheet.write_row(["UMI count normalization",normalization_method])
        worksheet.write_row(["Highly variable gene selection",hvg_method])
        worksheet.write_row(["Cell clustering","PCA and K-means clustering"])
        worksheet.write_row(["Differential expression analysis,SCDE"])
    
    workbook.close()
