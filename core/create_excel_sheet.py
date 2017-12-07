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
            #yield [ast.literal_eval(e) for e in line.strip('\n').split('\t')] ## Typecasting for right aligning numbers in excel

def write_excel_workbook(files_to_write,output_excel,catalog_number=None):
    ''' Write the give files as an excel workbook with different sheets
    :param list files_to_write: the list of tsv files to write
    :param str output_excel: the output excel file path
    :param str catalog_number: specify a catalog number, if applicable. The gene count excel sheet will be named so
    '''
    workbook = Workbook(output_excel)
    for infile in files_to_write:        
        if infile.find("gene") != -1 and catalog_number: ## Gene count file with a catalog number
            sheet_name = catalog_number
        else:
            sheet_name = os.path.basename(infile.rstrip("/"))
        worksheet = workbook.add_worksheet(sheet_name)
        i=0
        for row in file_reader(infile):
            worksheet.write_row(i,0,row)
            i+=1
    workbook.close()
