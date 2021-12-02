from biobookshelf.main import *

def Cigartuple_Convert( r_mappy = None, seq_mappy = None, cigartuple_pysam = None ) :
    """ # 2021-11-25 13:36:51 
    convert cigartuple from mappy to that of pysam and vice versa
    The output format will be fixed to list of tuples (pysam format)
    for mappy cigartuple -> pysam cigartuple, since the mappy does not include soft-clipping or hard clipping operations in its cigartuples, the hit record along with mappy sequence is required to return a complete cigartuples object (the clipped regions will be considered as soft-clipped regions)
    """
    # define interger representation of the CIGAR operations used in BAM files
    int_cigarop_S = 4
    int_cigarop_H = 5
    
    l_cigartuple = [ ]
    if r_mappy is not None and seq_mappy is not None :
        ''' mappy cigartuple -> pysam cigartuple '''
        r = r_mappy
        len_seq = len( seq_mappy )
        
        ''' handle soft clipping at the front (front in the '+' strand) '''
        if r.strand > 0 and r.q_st > 0 : # if sequence was aligned to '+' strand 
            l_cigartuple.append( ( int_cigarop_S, r.q_st ) )
        elif r.strand < 0 and r.q_en < len_seq :
            l_cigartuple.append( ( int_cigarop_S, len_seq - r.q_en ) )
            
        for int_n_base_pairs, int_cigarop in r.cigar :
            l_cigartuple.append( ( int_cigarop, int_n_base_pairs ) )
            
        ''' handle soft clipping at the end (end in the '+' strand) '''
        if r.strand > 0 and r.q_en < len_seq : # if sequence was aligned to '+' strand 
            l_cigartuple.append( ( int_cigarop_S, len_seq - r.q_en ) )
        elif r.strand < 0 and r.q_st > 0 :
            l_cigartuple.append( ( int_cigarop_S, r.q_st ) )
            
    elif cigartuple_pysam is not None :
        ''' pysam cigartuple -> mappy cigartuple '''
        for int_cigarop, int_n_base_pairs in cigartuple_pysam :
            l_cigartuple.append( ( int_n_base_pairs, int_cigarop ) )
    return l_cigartuple

def Generate_Base_and_Qual( r ) :
    """ # 2021-08-27 22:09:31 
    a generator function yielding the base and quality of every matched/mismatched positions of a given pysam read record 
    """
    pos_read, pos_ref = 0, 0 # current position in the extracted reference sequence and the current read
    int_oper_M, int_oper_I, int_oper_D, int_oper_N, int_oper_S, int_oper_H, int_oper_P, int_oper_equal, int_oper_X = list( range( 9 ) )
    for int_oper, n_bases in r.cigartuples :
        if int_oper == int_oper_M : 
            for _ in range( n_bases ) :
                str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
                yield r.reference_start + pos_ref, str_base_read, str_qual_read
                pos_ref += 1
                pos_read += 1
        elif int_oper == int_oper_N :
            pos_ref += n_bases
        elif int_oper == int_oper_S :
            pos_read += n_bases
        elif int_oper == int_oper_I :
            pos_read += n_bases
        elif int_oper == int_oper_D :
            pos_ref += n_bases
        elif int_oper == int_oper_H :
            pass
        elif int_oper == int_oper_P :
            pass
        elif int_oper == int_oper_equal :
            str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
            yield r.reference_start + pos_ref, str_base_read, str_qual_read
            pos_ref += n_bases
            pos_read += n_bases
        elif int_oper == int_oper_X :
            for _ in range( n_bases ) :
                str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
                yield r.reference_start + pos_ref, str_base_read, str_qual_read
                pos_ref += 1
                pos_read += 1
                
                
def Retrive_List_of_Mapped_Segments( cigartuples, pos_start, return_1_based_coordinate = False, flag_pos_start_0_based_coord = True, flag_return_splice_junction = False, flag_is_cigartuples_from_mappy = False ) :
    ''' # 2021-09-07 10:17:55 
    return l_seq and int_total_aligned_length for given cigartuples (returned by pysam cigartuples) and 'pos_start' (0-based coordinates, assuming pos_start is 0-based coordinate)
    'return_1_based_coordinate' : return 1-based coordinate, assuming 'pos_start' is 0-based coordinate (pysam returns 0-based coordinate)
    'flag_return_splice_junction' : additionally return list of splice junction tuples
    'flag_is_cigartuples_from_mappy' : if cigartuples are from mappy.Alignment, set this flag to True
    '''
    if return_1_based_coordinate and flag_pos_start_0_based_coord : # 0-based > 1-based
        start -= 1
    l_seg, start, int_aligned_length, int_total_aligned_length = list( ), pos_start, 0, 0
    for operation, length in cigartuples :
        if flag_is_cigartuples_from_mappy : # if flag_is_cigartuples_from_mappy, swap the two values
            temp = operation
            operation = length
            length = temp
        if operation in { 0, 2, 7, 8 } : # 'MD=X'
            int_aligned_length += length
        elif operation == 3 : # 'N' if splice junction appears, split the region and make a record
            l_seg.append( ( start, ( start + int_aligned_length - 1 ) if return_1_based_coordinate else ( start + int_aligned_length ) ) ) # set the end position
            start = start + int_aligned_length + length # set the next start position
            int_total_aligned_length += int_aligned_length # update total aligned length
            int_aligned_length = 0
    if int_aligned_length > 0 : 
        l_seg.append( ( start, ( start + int_aligned_length - 1 ) if return_1_based_coordinate else ( start + int_aligned_length ) ) )
        int_total_aligned_length += int_aligned_length
    if flag_return_splice_junction :
        # collect splice junction tuples from l_seg
        l_sj = [ ]
        for i in range( len( l_seg ) - 1 ) :
            l_sj.append( ( l_seg[ i ][ 1 ], l_seg[ i + 1 ][ 0 ] ) )
        return l_seg, int_total_aligned_length, l_sj
    else :
        return l_seg, int_total_aligned_length
    
def Alignment_Record_from_Mappy_to_String( hit, r_fastq ) :
    """ # 2021-10-27 15:19:48  
    For single-end read only
    
    'hit': mappy.Alignment object returned by mappy.Aligner.map function
    'r_fastq' : list of qname, seq, _, qual containing a single fastq record whose sequence has been queried using mappy
    """
    # set cigar operation values
    int_cigarop_S = 4
    int_cigarop_H = 5
    
    # parse fastq record
    header, seq, _, qual = r_fastq 
    qname = header[ 1 : ].split( ' ', 1 )[ 0 ] # retrieve qname
    
    flag_is_reverse_complemented = hit.strand != 1 
    int_flag = ( not hit.is_primary ) * 256 + ( flag_is_reverse_complemented ) * 16 # compose flag integer
    
    # retrieve seq and qual in the right orientation
    seq_in_sam_record = SEQ.Reverse_Complement( seq ) if flag_is_reverse_complemented else seq
    qual_in_sam_record = qual[ : : -1 ] if flag_is_reverse_complemented else qual
    
    # remove hard clipped portion of read (front)
    int_len_op, int_type_op = hit.cigar[ 0 ]
    if int_type_op == int_cigarop_H :
        seq_in_sam_record = seq_in_sam_record[ int_len_op : ]
        qual_in_sam_record = qual_in_sam_record[ int_len_op : ]
    # remove hard clipped portion of read (rear)
    int_len_op, int_type_op = hit.cigar[ -1 ]
    if int_type_op == int_cigarop_H :
        seq_in_sam_record = seq_in_sam_record[ : - int_len_op ]
        qual_in_sam_record = qual_in_sam_record[ : - int_len_op ]
        
    str_sam_record = '\t'.join( [ qname, str( int_flag ), hit.ctg, str( hit.r_st + 1 ), str( hit.mapq ), hit.cigar_str, '*', '0', '0', seq_in_sam_record, qual_in_sam_record ] ) # 0>1-based coordinates
    return str_sam_record

''' create a list of bookmarks for multithreading (simultaneous access of BAM files) '''
def Bookmark( dir_file_bam, n_threads, float_min_mapq_unique_mapped, int_n_bases_padding_around_interval = 10, verbose = False, l_dict_it = [ ], flag_ignore_record_with_flag_secondary_alignment = True, flag_ignore_record_with_flag_optical_or_pcr_duplicate = True ) :
    """  # 2021-11-10 23:35:29 
    iterate the given BAM file twice to create a bookmark to mark the start and end of each chunk for simultaneous access of the BAM file
    WARNING: there might be a little overlap of reads analyzed by each process since the chunk boundary is marked with ( qname, flag, refname, refstart ). If the two SAM records with the same ( qname, flag, refname, refstart ), the reads between the two 'identical' SAM records might be analyzed twice by two processes analyzing the two chunk sharing the boundary (however, it is unlikely that this will happens in a typical BAM file).
    
    'dir_file_bam' : indexed (with .bai file) BAM file to be create bookmarks for access by multiple processes
    'n_threads' : number of threads (or number of bookmarks to create) to access a single BAM file 
    'l_dict_it' : a list of dictionaries of interval trees (pip intervaltree package) containing intervals (genes and repeat annotations, etc.). When boundaries between chunk overlaps with the given intervals, it will be adjusted so that the boundary is at least 'int_n_bases_padding_around_interval' base pair away from the interval. 
        This functionality will allow analyzing all reads aligned to a single interval in the 'dict_it' to be analyzed in a single thread, thus abolishing the need to combine results from all threads to produce interval-level output, thus simplifying  downstream analysis.
    
    returns
    df_bookmark, l_int_n_sam_records_per_chunk
    """
    if isinstance( l_dict_it, dict ) : # if a single dict_it is given, wrap the dict_it inside a list
        l_dict_it = [ l_dict_it ]
    
    ''' retrieve sequence length of reference sequence from the BAM header '''
    with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
        dict_header = samfile.header.to_dict( )
    dict_seqname_to_len_seq = pd.DataFrame( dict_header[ 'SQ' ] ).set_index( 'SN' ).LN.to_dict( )
    
    
    """ Compose list of bookmarks for iterating input BAM file using multiple threads """
    l_col = [ 'id_refname_start', 'int_refstart_start', 'int_index_chunk_start', 'id_refname_end', 'int_refstart_end', 'int_index_chunk_end' ]
    # when number of threads > 1, retrieve bookmarks for multiprocessing
    if n_threads > 1 : 
        ''' retrieve the total number of SAM records '''
        # iterate the entire BAM file to retrieve the total number of SAM records
        int_n_total_sam_records = 0
        with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
            for r in samfile.fetch( ) :
                ''' filter alignments based on mapq and flag '''
                if r.mapq < float_min_mapq_unique_mapped : # skip read whose mapq is below 'float_min_mapq_unique_mapped'
                    continue
                # ignore records with flag indicating secondary alignment
                if flag_ignore_record_with_flag_secondary_alignment and ( r.flag & ( 1 << 8 ) ) :
                    continue
                # ignore records with flag indicating the alignment is optical pcr duplicates
                if flag_ignore_record_with_flag_optical_or_pcr_duplicate and ( r.flag & ( 1 << 10 ) ) :
                    continue
                ''' count the number of valid alignment records '''
                int_n_total_sam_records += 1
        
        """ retrieve the positions and the id_reads at the start and stop points ('bookmark') for multiprocessing """
        ''' retrieve bookmarks to mark the start and ends of contigs for each chunnk '''
        # calculate the number of reads for each chunk
        l_int_n_sam_records_per_chunk = INT.Split( int_n_total_sam_records, n_threads )

        l_bookmark = [ ] # empty bookmark
        flag_has_bookmark_been_initialized = False
        int_start_possible_boundary = None # start position of the possible boundary # will be set to non-None value when the desired number of records in the current chunk has been reached its target and the boundary is being set or being adjusted
        int_n_sam_record_count = 0
        int_index_chunk = 0
        id_refname_current = None # default reference name
        
        l_int_actual_n_sam_records_per_chunk = [ ] # actual number of sam records per chunk for the created bookmarks
        # iterate the entire BAM file to retrieve the total number of SAM records
        with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
            for r in samfile.fetch( ) :
                ''' filter alignments based on mapq and flag '''
                if r.mapq < float_min_mapq_unique_mapped : # skip read whose mapq is below 'float_min_mapq_unique_mapped'
                    continue
                # ignore records with flag indicating secondary alignment
                if flag_ignore_record_with_flag_secondary_alignment and ( r.flag & ( 1 << 8 ) ) :
                    continue
                # ignore records with flag indicating the alignment is optical pcr duplicates
                if flag_ignore_record_with_flag_optical_or_pcr_duplicate and ( r.flag & ( 1 << 10 ) ) :
                    continue
                    
                # retrieve SAM record info.
                refname = r.reference_name 
                refstart = r.reference_start
                refend = r.reference_end
                qname = r.qname
                flag = r.flag
                
                # initialize the bookmarks
                if not flag_has_bookmark_been_initialized :
                    l_bookmark.append( [ refname, refstart, 0 ] ) # initialize bookmark by assigning the first sam record to the index_chunk 0
                    id_refname_current = refname
                    flag_has_bookmark_been_initialized = True
                    continue
                
                int_n_sam_record_count += 1
                
                """ mark the change of contig by adding bookmarks """
                ''' if the next chromosome has been reached while for searching for next gene boundary, the boundary is considered passed, and the new chunk will be initiated.'''
                if id_refname_current != refname :
                    if verbose :
                        print( f'chromosome changed from {id_refname_current} to {refname} at chunk_{int_index_chunk}:{int_n_sam_record_count}' )
                    id_refname_current = refname
                    ''' handle the case when contig changes while the chunk boundary is being searched '''
                    if int_start_possible_boundary is not None : 
                        # retrieve stats about the current chunk and initializes the next chunk
                        int_start_possible_boundary = None
                        l_int_actual_n_sam_records_per_chunk.append( int_n_sam_record_count - 1 )
                        if verbose :
                            print( f'new chunk initiated (chunk {int_index_chunk + 1}) (chunk {int_index_chunk} has {int_n_sam_record_count} number of reads)' )
                        int_index_chunk += 1
                        int_n_sam_record_count = 1
                    bookmark = [ refname, 0, int_index_chunk ]
                    l_bookmark.append( bookmark ) # add bookmarks twice
                    l_bookmark.append( bookmark )

                """
                Search boundaries 
                """
                # initialize the next chunk if the 'int_n_sam_record_count' reached the assigned number of sam records for the current chunk (process)
                if int_n_sam_record_count > l_int_n_sam_records_per_chunk[ int_index_chunk ] :
                    ''' set default boundary (current position + 1) '''
                    if int_start_possible_boundary is None :
                        int_start_possible_boundary = refstart + 1
                    ''' iterate until reaching the possible chunk boundary '''
                    if int_start_possible_boundary is not None :
                        if refstart < int_start_possible_boundary : # if the next possible boundary has not been reached 
                            continue
                    
                    """
                    Possible boundary has been reached
                    """
                    ''' check whether the current position overlaps with intervals in the given list of dictionaries of interval trees, and set the next possible boundary based on the overlapped intervals '''
                    l_adjusted_boundary_based_on_overlapped_interval = [ ]
                    # search overlapped intervals
                    int_len_seq_current_contig = dict_seqname_to_len_seq[ refname ] # retrieve the length of the current contig
                    for dict_it in l_dict_it :
                        if refname in dict_it :
                            for start, end, anno_name in dict_it[ refname ][ refstart ] :
                                l_adjusted_boundary_based_on_overlapped_interval.append( min( end + int_n_bases_padding_around_interval, int_len_seq_current_contig - 1 ) ) # maximum boundary is the end position of the contig
                                if verbose :
                                    print( f"for chunk {int_index_chunk}, the possible boundary currently overlaps with {refname}:{start}-{end}|{anno_name}" )
                    # if there is any overlaps with the intervals, start adjusting boundaries
                    if len( l_adjusted_boundary_based_on_overlapped_interval ) != 0 :
                        int_start_possible_boundary_previous = int_start_possible_boundary
                        # set the next possible boundary
                        int_start_possible_boundary = max( l_adjusted_boundary_based_on_overlapped_interval ) # use the farthest interval ends when setting the next potential boundary 
                        if verbose :
                            print( f'for chunk {int_index_chunk}, the possible boundary will change from {refname}:{int_start_possible_boundary_previous} to {refname}:{int_start_possible_boundary} ({int_start_possible_boundary - refstart}bp shift) when int_n_sam_record_count = {int_n_sam_record_count}' )
                        continue
                        
                    """
                    Finalize boundary
                    """
                    ''' if the current position does not overlap with any intervals, set the current position as the boundary between chunks '''
                    int_start_possible_boundary = None # set 'int_start_possible_boundary' to indicate boundary is not being searched
                    l_int_actual_n_sam_records_per_chunk.append( int_n_sam_record_count - 1 )
                    if verbose :
                        print( f'new chunk initiated (chunk {int_index_chunk + 1}) (chunk {int_index_chunk} has {int_n_sam_record_count} number of reads)' )
                    # initialized the next chunk
                    int_n_sam_record_count = 1
                    int_index_chunk += 1
                    # add a bookmark of the current position to mark the start of the new chunk
                    bookmark = [ refname, refstart, int_index_chunk ]
                    l_bookmark.append( bookmark )
                    l_bookmark.append( bookmark )

        l_int_actual_n_sam_records_per_chunk.append( int_n_sam_record_count )
        l_bookmark.append( [ '', -1, int_index_chunk ] )
        l_bookmark = np.array( l_bookmark, dtype = object ).reshape( int( len( l_bookmark ) / 2 ), 6 ) # reshape bookmark so that each record contains bookmark of the start point and end points
        df_bookmark = pd.DataFrame( l_bookmark, columns = l_col )
    else :
        df_bookmark = pd.DataFrame( [ [ 'single_thread' ] * 6 ], columns = l_col )
    return df_bookmark, l_int_actual_n_sam_records_per_chunk