import regex

def endogenous_seq_match_clean(cigar,primer_len,read_is_reverse,num_flanking_bases=30,cigars_to_ignore=['N','D','H','P'], pattern=regex.compile('([0-9]+)([A-Z])')):
    ''' ## Check if query in gene
    '''

    expanded_cigar = []
    for num_bases,cigar_char in regex.findall(pattern,cigar):
        if cigar_char in cigars_to_ignore:
            continue
        else:
            expanded_cigar.extend(cigar_char*int(num_bases))

    loci1 = primer_len + num_flanking_bases + 1
    loci2 = primer_len

    if read_is_reverse:
        cigar_to_search = expanded_cigar[-(primer_len+num_flanking_bases):-(primer_len+1)]
    else:
        cigar_to_search = expanded_cigar[primer_len:(primer_len+num_flanking_bases+1)]

    return (cigar_to_search.count('M') >= 25)

def endogenous_seq_match(cigar,primer_len):
    '''
    '''

    j = -1
    query_counter = 0
    post_primer_matches = 0
    post_primer_bases = 0
    flag = 0

    for i in range(len(cigar)):
        if cigar[i].isalpha(): ## Store the index of cigar character
            num_bases = int(cigar[j+1:i])
            if cigar[i] in ['N','D','H','P']: ## These cigar characters do not consume the query
                j = i
                continue
            else:
                query_counter+=num_bases ## Increment query counter

            if query_counter > primer_len: ## Moved passed the primer part
                if flag == 0:## If this is the first time the primer region is being crossed
                    post_primer_bases += query_counter - primer_len
                else:
                    post_primer_bases += num_bases
                if cigar[i] == 'M':
                    if flag == 0:
                        post_primer_matches += query_counter - primer_len
                        flag = 1 ## set the flag
                    else:
                        post_primer_matches += num_bases

                ## Check if we have enough matches
                if post_primer_matches >= 20: ## Match for endogenous sequence
                    return True
                if post_primer_bases >= 30: ## Moved 30 bases past primer and still do not have an alignment of 20 bases
                    return False
            j = i

    return False

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

    read_sequence,read_is_reverse,read_len,read_chrom,read_pos,read_cigar,mt = read_tup
    if read_chrom == '*':
        return ('Unmapped','Unknown',0)
    if read_chrom not in primer_tree:
        return ('Unknown_Chrom','Unknown',0)

    if read_is_reverse: ## The primer is mapped to the last n bases of the read , as the reads in the bam are always on the +ve strand , hence we need to see if the end of the read still falls within the primer stop site from the design file.
        loci_to_search = read_pos + (read_len-2)
    else:
        loci_to_search = read_pos

    res = primer_tree[read_chrom].search(loci_to_search-1,loci_to_search+1) ## allowing a shift in primer start loci by 1 base pair , account for soft clip in the begining 
    if res:
        if len(res) > 1:
            raise Exception("error in primer finding, read spans more than 2 primer loci") ## take multiple hits
        else:
            result = res.pop()
            pattern,primer = result.data
            if regex.match(pattern,read_sequence): ## Check if the primer has approximate match to the read sequence
                ## Check for endogenous sequence
                if endogenous_seq_match_clean(read_cigar,len(primer),read_is_reverse):
                    return (primer , mt, 1)
                else:
                    return (primer, mt, 0)
            else:
                return ('Unknown_Regex','Unknown',0)
    else:
        return('Unknown_Loci','Unknown',0)
    
    #return (primer,mt,0) ## mispriming : either the loci was not on any valid primer target or the regex match failed
