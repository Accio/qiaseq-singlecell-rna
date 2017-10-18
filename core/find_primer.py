import regex

def endogenous_seq_match(cigar,primer_len,read_is_reverse,num_flanking_bases=30,cigars_to_ignore=['N','D','H','P'], pattern=regex.compile('([0-9]+)([A-Z])')):
    ''' Function to check whether the sequence following the primer has enough matching alignment
    :param str cigar: The cigar string , e.g. 10S80M1D200N20M
    :param int primer_len: The length of the primer
    :param bool read_is_reverse: Whether the read is reverse complemented
    :param int num_flanking_bases: The number of flanking bases after the primer to check the cigar for matches (default)
    :param list cigars_to_ignore: the list of cigar characters to ignore when walking along the cigarstring, these should be the characters which do not consume the query sequence (default)
    :param regex object patter: the regex pattern for the cigar string (default)    
    '''
    ## Expand the cigar string and store it in a list
    expanded_cigar = []
    for num_bases,cigar_char in regex.findall(pattern,cigar):
        if cigar_char in cigars_to_ignore:
            continue
        else:
            expanded_cigar.extend(cigar_char*int(num_bases))
    ## Choose appropriate coordinates to check for the cigar matches
    loci1 = primer_len + num_flanking_bases + 1
    loci2 = primer_len
    ## Match the cigar based on strandedness
    if read_is_reverse:
        cigar_to_search = expanded_cigar[-(primer_len+num_flanking_bases):-(primer_len+1)]
    else:
        cigar_to_search = expanded_cigar[primer_len:(primer_len+num_flanking_bases+1)]

    return (cigar_to_search.count('M') >= 25)

def find_primer(primer_tree,read_tup):
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
    read_name,read_sequence,read_is_reverse,read_len,read_chrom,read_pos,read_cigar,mt,nh = read_tup
    if read_chrom == '*':
        return ('Unmapped',mt,0,0)
    if read_chrom not in primer_tree:
        return ('Unknown_Chrom',mt,0,nh)

    if read_is_reverse: ## The primer is mapped to the last n bases of the read , as the reads in the bam are always on the +ve strand , hence we need to see if the end of the read still falls within the primer stop site from the design file.
        loci_to_search = read_pos + (read_len-1)
        res = primer_tree[read_chrom].search(loci_to_search,loci_to_search-3) ## allowing some offset in primer start loci 
    else:
        loci_to_search = read_pos
        res = primer_tree[read_chrom].search(loci_to_search,loci_to_search+3)

    if res:
        for i in range(len(res)):
            result = res.pop()
            pattern,primer = result.data
            if regex.match(pattern,read_sequence): ## Check if the primer has approximate match to the read sequence
                ## Check for endogenous sequence
                if endogenous_seq_match(read_cigar,len(primer),read_is_reverse):
                    return (primer , mt, 1,nh)
                else:
                    return (primer, mt, 0,nh)            
        return ('Unknown_Regex',mt,0,nh)
    else:
        return('Unknown_Loci',mt,0,nh)
