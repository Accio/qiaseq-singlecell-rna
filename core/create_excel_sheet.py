import os
from xlsxwriter.workbook import Workbook
import ast

def file_reader(infile):
    ''' Reads a file and yields the file line
    :param str infile: a tsv file    
    yields the rows of the file as a list
    '''
    with open(infile,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')
            type_casted_row = []
            for e in contents:
                try: ## Try int
                    val = int(e)
                except ValueError: ## Try float
                    try:
                        val = float(e)                    
                    except ValueError: ## Keep as string
                        val = e
                type_casted_row.append(val)
            yield type_casted_row

def write_excel_workbook(files_to_write,output_excel,catalog_number=None,species=None):
    ''' Write the give files as an excel workbook with different sheets
    :param list files_to_write: the list of tsv files to write
    :param str output_excel: the output excel file path
    :param str catalog_number: specify a catalog number, if applicable. The gene count excel sheet will be named so
    :param str catalog_number: specify a species name, if applicable. The gene count excel sheet will be have this
    '''
    workbook = Workbook(output_excel)
    for infile in files_to_write:
        if infile.find("gene") != -1 and catalog_number: ## Gene count file with a catalog number
            sheet_name = "umis.genes."+catalog_number
        elif infile.find("gene") != -1 and species: ## Gene count file with a species (i.e. polyA)
            sheet_name = "umis.genes.polyA-"+species.lower()
        elif infile.find("primer") != -1: ## Primer count file
            sheet_name = "umis.primers."+catalog_number
        elif infile.find(".metrics.by_sample_index") != -1: ## Sample Index metrics
            sheet_name = "metrics.by_sample_index"
        elif infile.find(".metrics.by_cell_index") != -1: ## Cell Index metrics
            sheet_name = "metrics.by_cell_index"
        else:
            raise Exception("Invalid file name encountered !")
        worksheet = workbook.add_worksheet(sheet_name)
        i=0
        for row in file_reader(infile):
            worksheet.write_row(i,0,row)
            i+=1
    workbook.close()
