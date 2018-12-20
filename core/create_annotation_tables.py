import sqlite3
import os
import io
import gzip
from intervaltree import IntervalTree
from collections import defaultdict


def open_by_magic(filename):
    '''
    Adapted from : http://stackoverflow.com/questions/18367511/how-do-i-automatically-handle-decompression-when-reading-a-file-in-python
    with modifications
    Uses the initial bytes of a file to detect the file compression.

    :param str filename: path to the input file
    :return: the appropriate file handle for reading
    :rtype: file object
    '''
    ## Add more magic strs here for various compressions
    magic_dict = {"\x1f\x8b\x08":gzip.open}
    max_len = max(len(x) for x in magic_dict)
    with open(filename) as f:
        file_start = f.read(max_len)
        for magic,fn in magic_dict.items():
            if file_start.startswith(magic):
                return io.BufferedReader(fn(filename))
            return open(filename,'r') ## Otherwise just a regular file

def convert_strand(strand):
    ''' Convert strand to QIAGEN format
    :param str strand: the strand , either + or -
    :return '1' or '-1'
    :rtype str
    '''
    if strand == "+":
        return '1'
    elif strand == '-':
        return '-1'
    else:
        raise Exception("Invalid Strand information !")

def create_gene_hash(annotation_gtf,ercc_bed):
    ''' Create a Hash table with annotation information for Genes

    :param str: annotation_gtf: path to the genic annotation gtf file 
    :param str: ercc_bed: path to the bed file with ERCC information

    :return A dictionary with annotation information
    :rtype dict
    '''
    gene_info = defaultdict(list)
    ## Parse the gencode annotation file and store in the dictionary
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            if line[0] == "#":
                continue
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                start = contents[3]
                end = contents[4]
                strand = convert_strand(contents[6])                    
                info = contents[-1]
                gene = info.split(';')[3].split()[1].strip('\"')
                gene_type = info.split(';')[1].split()[1].strip('\"')
                cols = [chrom,start,end,strand,gene,gene_type]
                gene_info[gene] = cols
                
    with open(ercc_bed,'r') as IN:
        for line in IN:
            chrom,start,end,seq,strand,ercc = line.strip("\n").split("\t")
            cols = [chrom,start,end,strand,chrom,ercc]
            gene_info[chrom] = cols

    return gene_info

def create_gene_tree(annotation_gtf,ercc_bed,species,merge_coordinates=False):
    '''
    :param str annotation_gtf : a gtf file for identifying genic regions
    :param str ercc_bed: a bed file for storing information about ERCC regions
    :param str species: species name (qiagen's internal alias)
    
    :return An interval tree with annotation gene and ercc annotation information
    :rtype object : IntervalTree data structure
    '''
    gene_tree = defaultdict(lambda:defaultdict(IntervalTree))
    genes = defaultdict(list)
    valid_chromosomes = ["chr"+str(i) for i in range(0,23)]    
    valid_chromosomes.extend(["chrX","chrY","chrM","chrMT"])
    well_annotated = ['human','mouse']
    
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            if line[0]=='#':
                continue
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                if species.lower() in well_annotated and (chrom not in valid_chromosomes and not chrom.startswith('ERCC')): ## Skip contigs
                    continue
                start = int(contents[3])
                end = int(contents[4])
                strand = convert_strand(contents[6])
                info = contents[-1]
                gene,ensemb_id,gene_type = (None,None,None)
                for element in info.split(';'):
                    if len(element) == 0:
                        continue                    
                    e1,e2 = element.split()
                    if e1 == "gene_name":
                        gene = e2.strip('\"')
                    elif e1 in ["gene_type","gene_biotype"]:
                        gene_type = e2.strip('\"')
                    elif e1 == "gene_id":
                        ensembl_id = e2.strip('\"')               
                        
                if ensembl_id == None:
                    raise Exception(
                        "Failed Parsing annotation file :{annotation}".format(
                            annotation=annotation_gtf))
                if gene == None: ## Ensembl rat gtf did not have gene name in certain cases
                    gene = ensembl_id

                if merge_coordinates: ## Create a coordinate set which is merged to include the largest interval possible
                    if gene in genes: ## Gene has been seen before
                        if genes[gene][2] != chrom: ## Different chromosome
                            genes[gene] = (start,end,chrom,strand,gene,ensembl_id)
                        else:
                            if start <= genes[gene][0]: ## Gene has unique start bases to add
                                if end >= genes[gene][1]: ## Gene has unique end bases to add
                                    genes[gene] = (start,end,chrom,strand,gene,ensembl_id)
                                else:
                                    if end < genes[gene][0]: ## Check if this interval is dijoint from the previous one
                                        genes[gene].append((start,end,chrom,strand,gene,ensembl_id)) ## Add a new interval
                                    else:
                                        genes[gene] = (start,genes[gene][1],chrom,strand,gene,ensembl_id) ## Merge intervals

                            else: ## Start is already spanned , check end
                                if end >= genes[gene][1]: ## Update end base position
                                    if start > genes[gene][1]: ## Start is greater than previously encountered gene's end
                                        genes[gene].append((start,end,chrom,strand,gene,ensembl_id)) ## Add a new interval
                                    else:
                                        genes[gene] = (genes[gene][0],end,chrom,strand,gene,ensembl_id) ## Merge intervals
                                else: ## No need to update anything
                                    continue
                    else:
                        genes[gene] = (start,end,chrom,strand,gene,ensembl_id)

                else: ## Store all intervals without merging
                    genes[gene].append((start,end,chrom,strand,gene,ensembl_id))
    
    ## Build a gene tree to store gene info
    for gene in genes:
        for info in genes[gene]:
            start,end,chrom,strand,gene,ensembl_id = info
            assert strand in ['-1','1'],"Incorrect strand !"
            if strand == "1":
                five_prime = start
                three_prime = end
            else:
                five_prime = end
                three_prime = start
            new_info = (ensembl_id,gene,strand,chrom,five_prime,three_prime)
            gene_tree[chrom][strand].addi(start,end+1,new_info)

    print "Interval tree created with {ngenes} genes".format(ngenes=len(genes))
    return gene_tree
