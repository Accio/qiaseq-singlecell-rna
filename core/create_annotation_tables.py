import sqlite3
import os
from intervaltree import IntervalTree
from collections import defaultdict
from extract_multiplex_region import open_by_magic

def create_gene_anno(annotation_db_file,annotation_gtf):
    ''' Create a table to store genic annotations
    '''
    ## Remove the file if it is already present
    os.system('rm -f {}'.format(annotation_db_file))
    conn = sqlite3.connect(annotation_db_file)
    cur = conn.cursor()
    create_table_sql = (
        '''CREATE TABLE gene_annotation '''
        '''(chrom text,start int,stop int,strand text'''
        ''',info text,gene text,gene_type text)'''
    )
    cur.execute(create_table_sql)
    ## Parse the gencode annotation file and insert to db
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            if line[0] == "#": ## Skip header
                continue
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                start = contents[3]
                end = contents[4]
                strand = contents[6]
                info = contents[-1]
                gene = info.split(';')[3].split()[1].strip('\"')
                gene_type = info.split(';')[1].split()[1].strip('\"')
                cols = [chrom,start,end,strand,info,gene,gene_type]
                cur.execute('''INSERT INTO gene_annotation VALUES (?,?,?,?,?,?,?)''',cols)
                conn.commit()
    conn.close()

def create_gene_hash(annotation_gtf):
    '''
    '''
    gene_info = defaultdict(list)
    ## Parse the gencode annotation file and insert to db
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            if line[0] == "#":
                continue
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                start = contents[3]
                end = contents[4]
                strand = contents[6]
                info = contents[-1]
                gene = info.split(';')[3].split()[1].strip('\"')
                gene_type = info.split(';')[1].split()[1].strip('\"')
                cols = [chrom,start,end,strand,gene,gene_type]
                #if gene in gene_info:
                    #print 'duplicate gene: {}'.format(gene)
                gene_info[gene] = cols
    return gene_info

def create_gene_tree(annotation_gtf,merge_coordinates=False):
    '''
    :param str annotation_gtf : a gtf file for identifying genic regions
    '''
    gene_tree = defaultdict(lambda:defaultdict(IntervalTree))
    genes = defaultdict(list)
    with open_by_magic(annotation_gtf) as IN:
        for line in IN:
            if line[0]=='#':
                continue
            contents = line.strip('\n').split('\t')
            if contents[2] == 'gene':
                chrom = contents[0]
                start = int(contents[3])
                end = int(contents[4])
                strand = contents[6]
                info = contents[-1]
                gene = info.split(';')[3].split()[1].strip('\"')
                gene_type = info.split(';')[1].split()[1].strip('\"')

                if gene == None or gene_type == None:
                    raise Exception(
                        "Failed Parsing annotation file :{annotation}".format(
                            annotation=annotation_gtf))

                if merge_coordinates: ## Create a coordinate set which is merged to include the largest interval possible
                    if gene in genes: ## Gene has been seen before
                        if genes[gene][2] != chrom: ## Different chromosome
                            genes[gene] = (start,end,chrom,strand,gene,gene_type)
                        else:
                            if start <= genes[gene][0]: ## Gene has unique start bases to add
                                if end >= genes[gene][1]: ## Gene has unique end bases to add
                                    genes[gene] = (start,end,chrom,strand,gene,gene_type)
                                else:
                                    if end < genes[gene][0]: ## Check if this interval is dijoint from the previous one
                                        genes[gene].append((start,end,chrom,strand,gene,gene_type)) ## Add a new interval
                                    else:
                                        genes[gene] = (start,genes[gene][1],chrom,strand,gene,gene_type) ## Merge intervals

                            else: ## Start is already spanned , check end
                                if end >= genes[gene][1]: ## Update end base position
                                    if start > genes[gene][1]: ## Start is greater than previously encountered gene's end
                                        genes[gene].append((start,end,chrom,strand,gene,gene_type)) ## Add a new interval
                                    else:
                                        genes[gene] = (genes[gene][0],end,chrom,strand,gene,gene_type) ## Merge intervals
                                else: ## No need to update anything
                                    continue
                    else:
                        genes[gene] = (start,end,chrom,strand,gene,gene_type)

                else: ## Store all intervals without merging
                    genes[gene].append((start,end,chrom,strand,gene,gene_type))

    ## Build a gene tree to store gene info
    for gene in genes:
        for info in genes[gene]:
            start,end,chrom,strand,gene,gene_type = info
            gene_tree[chrom][strand].addi(start,end+1,info)

    print "Interval tree created with {ngenes}".format(ngenes=len(genes))
    return gene_tree
