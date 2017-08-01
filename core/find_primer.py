import regex

def find_primer(primer_dict,patterns,read_tup):
    '''
    Find whether a read matches one of the SPE primers used for the
    sequencing experiment

    :param dict of dict of lists primer_dict: Contains primer info
    :param a dict of dict  of lists patterns: Contains the compiled regular expression
                                              for the primer
    :param tuple read_tup: (read_sequence,chromosome,MT)
    :returns: a tuple containing the primer and mt info and whether it was a match
    :rtype: tuple
    '''

    read_sequence,read_chrom,mt = read_tup
    if read_chrom not in primer_dict:
        print "Invalid Chromosome for primer : %s"%(read_chrom)
        return ('Unknown','Unknown',0)

    for primer in primer_dict[read_chrom]:
        if regex.match(patterns[read_chrom][primer],read_sequence):
            return (primer , mt, 1)

    ## Fall back search on other chromosomes, making the search on other chromosomes to be
    ## strict , i.e. exact match
    for chrom in primer_dict:
        if chrom == read_chrom:
            continue
        for primer in primer_dict[chrom]:
            #if regex.match(patterns[chrom][primer],read_sequence):
            strand = primer_dict[chrom][primer][-2]
            if strand == '-':
                primer_seq = primer_dict[chrom][primer][-3]
            else:
                primer_seq = primer
            if primer_seq in read_sequence:
                return (primer, mt, 1)
            elif read_sequence in primer_seq:
                return(primer,mt,1)

    return (primer,mt,0)

