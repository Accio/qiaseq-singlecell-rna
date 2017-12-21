import regex
import logging
import RemoteException

def overlap(x1,x2,y1,y2):
    ''' Compute overlap between two coordinates
    '''
    if y1 > y2: ## Swap for genes on -ve strand
        temp = y2
        y2 = y1
        y1 = temp        
    return max(0,min(x2,y2) - max(x1,y1))

def return_read_end_pos(read_pos,cigar,flag=False,pattern=regex.compile('([0-9]+)([A-Z])'),cigars_to_ignore=['I','S','H','P']):
    ''' Return the end of the read
    '''
    expanded_cigar = []
    bases=0
    if flag:
        bases_to_ignore = ['S','I','H','P','N']
    for num_bases,cigar_char in regex.findall(pattern,cigar):
        if cigar_char in cigars_to_ignore:
            continue
    else:
        bases+=int(num_bases)
    return read_pos+bases

@RemoteException.showError
def find_gene(gene_tree,read_tup):
    ''' Annotate the given read with a gene

    :param interval treee gene_tree: efficient interval tree data structure storing coordinates and gene information
    :param tuple read_tup: a tuple of read information
    :return a tuple containing gene and mt info
    :rtype tuple
    '''

    overlap_threshold = 0
    logger = logging.getLogger("count_umis")
    read_id,read_sequence,read_is_reverse,read_len,read_chrom,read_pos,read_cigar,mt,nh = read_tup
    if read_chrom == "*":
        logger.info("{read_id}: Unmapped".format(read_id=read_id))
        return ('Unmapped',mt,0,nh)
    if 'ERCC' in read_chrom:
        read_end = return_read_end_pos(read_pos,read_cigar)
        result = gene_tree[read_chrom]["1"].search(read_pos,read_end)
        if result:
            return (next(iter(result)).data,mt,0,nh)            
            logger.info("{read_id}: Mapped to {ercc}".format(read_id=read_id,ercc=read_chrom))
        else:
            return ('Unknown_Chrom',mt,0,nh)
            logger.info("{read_id}: Mapped to {ercc}, but was off-loci in gene-tree".format(read_id=read_id,ercc=read_chrom))

    if read_chrom not in gene_tree:
        logger.info("{read_id}: Chromosome {chrom} was not present in annotation gene interval".format(read_id=read_id,chrom=read_chrom))
        return ('Unknown',mt,0,nh)

    ## Search the interval tree
    read_end = return_read_end_pos(read_pos,read_cigar)
    if read_is_reverse:
        res = gene_tree[read_chrom]['-1'].search(read_pos,read_end)
    else:
        res = gene_tree[read_chrom]['1'].search(read_pos,read_end)

    if res: ## If the search was successful
        num_hits = len(res)
        logger.info("{}".format(num_hits))
        if num_hits > 1:
            logger.info("{read_id}: intersected with the genes : \n{hits}".format(read_id=read_id,hits=res))
            ## Choose the closest 3' location gene
            prev = None
            for result in res:
                three_prime = result.data[4]
                five_prime = result.data[5]
                gene = result.data[1]
                if not prev: ## Checking first hit
                    ## Check overlap                    
                    prev_o = float(overlap(read_pos,read_end,five_prime,three_prime))/read_len
                    if prev_o > overlap_threshold:
                        prev = result.data
                        logger.info("{read_id}: Picked {gene} as default".format(read_id=read_id,gene=gene))
                    else:
                        logger.info("{read_id}: {gene} had {overlap} overlap, failed overlap criteria".format(read_id=read_id,overlap=prev_o,gene=gene))
                        
                else: ## Check other hits
                    o = float(overlap(read_pos,read_end,five_prime,three_prime))/read_len
                    if o < overlap_threshold:
                        continue
                    
                    three_prime_prev = prev[4]                        
                    three_prime_read = read_pos                   
                    diff_three_prime_prev = abs(three_prime_prev - three_prime_read)
                    diff_three_prime_current = abs(three_prime - three_prime_read)

                    if diff_three_prime_prev == diff_three_prime_current:
                        if o < prev_o: ## Look at overlaps
                            logger.info("{read_id}: Picked {gene1}: 3_prime_diff={diff1} over {gene2}: 3_prime_diff={diff2} because of greater overlap".format(read_id=read_id,gene1=gene,diff1=diff_three_prime_current,gene2=prev[1],diff2=diff_three_prime_prev))
                            prev_o = o
                            prev = result.data
                    elif diff_three_prime_prev > diff_three_prime_current:
                        logger.info("{read_id}: Picked {gene1}: 3_prime_diff={diff1} over {gene2}: 3_prime_diff={diff2}".format(read_id=read_id,gene1=gene,diff1=diff_three_prime_current,gene2=prev[1],diff2=diff_three_prime_prev))
                        prev = result.data
                        prev_o = o
            if prev:
                logger.info("{read_id}: Picked {gene1}".format(read_id=read_id,gene1=prev[1]))
                return (prev,mt,1,nh)
            else:
                logger.info("{read_id}: No genes matched overlap criteria".format(read_id=read_id))
                return ('Unknown',mt,0,nh)
        else:
            result = res.pop()
            logger.info("{read_id}: intersected with {gene} only".format(read_id=read_id,gene=result.data[1]))
            o = float(overlap(read_pos,read_end,result.data[5],result.data[4]))/read_len
            if o < overlap_threshold:
                logger.info("{read_id}: {gene} failed overlap criteria".format(read_id=read_id,gene=result.data[1]))
                return ('Unknown',mt,0,nh)
            else:
                return (result.data,mt,1,nh)
    else: ## Could not find loci in gene tree
        logger.info("{read_id} was not found in the annotation gene interval tree".format(read_id=read_id))
        return ('Unknown',mt,0,nh)
