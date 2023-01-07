from biobookshelf.main import *
import biobookshelf.SEQ as SEQ

def combine_bam( path_file_bam_output, * l_bam_file, int_max_num_files = 900, flag_remove_file = False, int_num_threads = 12 ) :
    """ # 2022-07-21 14:49:08 CML & HSA 
    abstraction: hide the detail implementation to users.

    'path_file_bam_output' : path to output file. if the output file already exists, this throw an error
    'l_bam_file' : list of bam files that will be merged
    'int_max_num_files' : maximum number of files that can be opened in a process (typicaly 1024 in linux)
    'flag_remove_file' : if True, input files will be removed.
    'int_num_threads' : the number of threads to merge input bam files. one thread will be used to distribute works (so int_num_threads - 1 threads will be actually used for sorting)
    """
    int_num_files = len( l_bam_file )
    int_num_of_recursive_merging = int( np.ceil( math.log( int_num_files, int_max_num_files ) ) ) # retrieve the minmum number of recursive merging

    path_folder_output = path_file_bam_output.rsplit( '/', 1 )[ 0 ] + '/' # retrieve output folder where the output file will reside
    path_folder_temp = f'{path_folder_output}temp_{UUID( )}/'
    
    # create folders
    for path_folder in [ path_folder_output, path_folder_temp ] :
        os.makedirs( path_folder, exist_ok = True ) # create temp folder

    # merge files until all files are merged in to a single file
    while len( l_bam_file ) > 1 : # while loop is active when there is files to be merged
        int_num_files = len( l_bam_file )  # retreive number of bam files
        l_file_size = list( os.stat( path_file ).st_size for path_file in l_bam_file ) # retrieve file size of input bam files

        str_uuid_current_merge_step = UUID( ) # identifier for the current merge step

        # retrieve number of output bam files
        if int_num_files < int_max_num_files :
            int_num_output_bam_file = 1
        else :
            int_num_output_bam_file = int( np.ceil( int_num_files / int_max_num_files ) ) 

        def batch_generator(  ) :
            l_batch = LIST_Split( l_bam_file, n_split = int_num_output_bam_file, flag_contiguous_chunk = True, arr_weight_for_load_balancing = l_file_size )
            for batch in l_batch :
                yield batch
        def process_batch( batch, pipe_sender = None ) :
            l_file_bam_to_be_merged_into_a_single_bam_file = batch 
            path_file_bam_combined_for_a_batch = f"{path_folder_temp}{str_uuid_current_merge_step}.{UUID( )}.bam"
            print( len( l_file_bam_to_be_merged_into_a_single_bam_file ), 'files will be merged into', path_file_bam_combined_for_a_batch )
            pysam.merge( '--no-PG', path_file_bam_combined_for_a_batch, * l_file_bam_to_be_merged_into_a_single_bam_file )
            pysam.index( path_file_bam_combined_for_a_batch )
            pipe_sender.send( 'completed' ) # tell the distributor that the work is done! 

        Multiprocessing_Batch( batch_generator(  ), process_batch, int_num_threads = max( 2, int_num_threads ) )

        # remove temporary or input files
        if flag_remove_file :
            for path_file in l_bam_file :
                os.remove( path_file )
                os.remove( path_file + '.bai' ) # delete index files, too
                pass
        flag_remove_file = True # remove temporary files

        l_bam_file = glob.glob( f"{path_folder_temp}{str_uuid_current_merge_step}.*.bam" ) # retrieve input bam files for the next round

    # rename combined output file
    os.rename( l_bam_file[ 0 ], path_file_bam_output )
    os.rename( l_bam_file[ 0 ] + '.bai', path_file_bam_output + '.bai' )

    # delete temporary folder
    shutil.rmtree( path_folder_temp )

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
    ''' # 2022-01-28 01:51:16 
    return l_seq and int_total_aligned_length for given cigartuples (returned by pysam cigartuples) and 'pos_start' (0-based coordinates, assuming pos_start is 0-based coordinate)
    'return_1_based_coordinate' : return 1-based coordinate, assuming 'pos_start' is 0-based coordinate (pysam returns 0-based coordinate)
    'flag_return_splice_junction' : additionally return list of splice junction tuples
    'flag_is_cigartuples_from_mappy' : if cigartuples are from mappy.Alignment, set this flag to True
    '''
    l_seg, start, int_aligned_length, int_total_aligned_length = list( ), pos_start, 0, 0
    if return_1_based_coordinate and flag_pos_start_0_based_coord : # 0-based > 1-based
        start -= 1
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

"""
Build a list of bookmarks for multithreading (simultaneous access of BAM files) 
"""
def Bookmark( path_file_bam, n_threads, int_min_mapq_unique_mapped, int_n_bases_padding_around_interval = 10, verbose = False, l_dict_it = [ ], flag_ignore_record_with_flag_secondary_alignment = True, flag_ignore_record_with_flag_optical_or_pcr_duplicate = True, flag_genome_without_chr_prefix_was_used_for_l_dict_it = True ) :
    """  # 2021-11-10 23:35:29 
    iterate the given BAM file twice to create a bookmark to mark the start and end of each chunk for simultaneous access of the BAM file
    WARNING: there might be a little overlap of reads analyzed by each process since the chunk boundary is marked with ( qname, flag, refname, refstart ). If the two SAM records with the same ( qname, flag, refname, refstart ), the reads between the two 'identical' SAM records might be analyzed twice by two processes analyzing the two chunk sharing the boundary (however, it is unlikely that this will happens in a typical BAM file).

    'path_file_bam' : indexed (with .bai file) BAM file to be create bookmarks for access by multiple processes
    'n_threads' : number of threads (or number of bookmarks to create) to access a single BAM file 
    'l_dict_it' : a list of dictionaries of interval trees (pip intervaltree package) containing intervals (genes and repeat annotations, etc.). When boundaries between chunk overlaps with the given intervals, it will be adjusted so that the boundary is at least 'int_n_bases_padding_around_interval' base pair away from the interval. 
        This functionality will allow analyzing all reads aligned to a single interval in the 'dict_it' to be analyzed in a single thread, thus abolishing the need to combine results from all threads to produce interval-level output, thus simplifying  downstream analysis. if chromosome name contains 'chr' prefix, please set 'flag_genome_without_chr_prefix_was_used_for_l_dict_it' to False
    'flag_genome_without_chr_prefix_was_used_for_l_dict_it' : (default = True) Please set this flag to False if the chromosome names of the given 'l_dict_it' annotations have a 'chr' prefix.

    returns
    df_bookmark, l_int_n_sam_records_per_chunk
    """
    if isinstance( l_dict_it, dict ) : # if a single dict_it is given, wrap the dict_it inside a list
        l_dict_it = [ l_dict_it ]

    ''' retrieve sequence length of reference sequence from the BAM header '''
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        dict_header = samfile.header.to_dict( )
    dict_seqname_to_len_seq = pd.DataFrame( dict_header[ 'SQ' ] ).set_index( 'SN' ).LN.to_dict( )

    """ define internal functions """
    def __chromosome_name_remove_chr__( str_name_chrom ) :
        ''' remove 'chr' prefix from the chromosome name '''
        if 'chr' == str_name_chrom[ : 3 ] :
            return str_name_chrom[ 3 : ]
        else :
            return str_name_chrom

    """ Compose list of bookmarks for iterating input BAM file using multiple threads """
    l_col = [ 'id_refname_start', 'int_refstart_start', 'int_index_chunk_start', 'id_refname_end', 'int_refstart_end', 'int_index_chunk_end' ]
    # when number of threads > 1, retrieve bookmarks for multiprocessing
    if n_threads > 1 : 
        ''' retrieve the total number of SAM records '''
        # iterate the entire BAM file to retrieve the total number of SAM records
        int_n_total_sam_records = 0
        with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
            for r in samfile.fetch( ) :
                ''' filter alignments based on mapq and flag '''
                if r.mapq < int_min_mapq_unique_mapped : # skip read whose mapq is below 'int_min_mapq_unique_mapped'
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
        with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
            for r in samfile.fetch( ) :
                ''' filter alignments based on mapq and flag '''
                if r.mapq < int_min_mapq_unique_mapped : # skip read whose mapq is below 'int_min_mapq_unique_mapped'
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
                    refname_for_l_dict_it = __chromosome_name_remove_chr__( refname ) if flag_genome_without_chr_prefix_was_used_for_l_dict_it else refname # handle if BAM file's genome chromosome name convention and 'l_dict_it' annotations' chromosome name convention is different # retrieve refname without 'chr' prefix if 'flag_genome_without_chr_prefix_was_used_for_l_dict_it' is True

                    for dict_it in l_dict_it :
                        if refname_for_l_dict_it in dict_it :
                            for start, end, anno_name in dict_it[ refname_for_l_dict_it ][ refstart ] :
                                l_adjusted_boundary_based_on_overlapped_interval.append( min( end + int_n_bases_padding_around_interval, int_len_seq_current_contig - 1 ) ) # maximum boundary is the end position of the contig
                                if verbose :
                                    print( f"for chunk {int_index_chunk}, the possible boundary currently overlaps with {refname_for_l_dict_it}:{start}-{end}|{anno_name}" )
                    # if there is any overlaps with the intervals, start adjusting boundaries
                    if len( l_adjusted_boundary_based_on_overlapped_interval ) != 0 :
                        int_start_possible_boundary_previous = int_start_possible_boundary
                        # set the next possible boundary
                        int_start_possible_boundary = max( l_adjusted_boundary_based_on_overlapped_interval ) # use the farthest interval ends when setting the next potential boundary 
                        if verbose :
                            print( f'for chunk {int_index_chunk}, the possible boundary will change from {refname_for_l_dict_it}:{int_start_possible_boundary_previous} to {refname_for_l_dict_it}:{int_start_possible_boundary} ({int_start_possible_boundary - refstart}bp shift) when int_n_sam_record_count = {int_n_sam_record_count}' )
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

def Bookmark_Generator( path_file_bam, int_min_mapq_unique_mapped, int_num_sam_records_for_each_chunk = 1000000, int_n_bases_padding_around_interval = 10, verbose = False, l_dict_it = [ ], flag_ignore_record_with_flag_secondary_alignment = True, flag_ignore_record_with_flag_optical_or_pcr_duplicate = True, flag_genome_without_chr_prefix_was_used_for_l_dict_it = True ) :
    """  # 2022-01-17 21:50:42 
    iterate the given BAM file twice to create a bookmark to mark the start and end of each chunk for simultaneous access of the BAM file.
    yield a single bookmark at a time
    WARNING: there might be a little overlap of reads analyzed by each process since the chunk boundary is marked with ( qname, flag, refname, refstart ). If the two SAM records with the same ( qname, flag, refname, refstart ), the reads between the two 'identical' SAM records might be analyzed twice by two processes analyzing the two chunk sharing the boundary (however, it is unlikely that this will happens in a typical BAM file).

    'path_file_bam' : indexed (with .bai file) BAM file to be create bookmarks for access by multiple processes
    'int_num_sam_records_for_each_chunk' : a guidance to the number of sam records contained in each chunk 
    'l_dict_it' : a list of dictionaries of interval trees (pip intervaltree package) containing intervals (genes and repeat annotations, etc.). When boundaries between chunk overlaps with the given intervals, it will be adjusted so that the boundary is at least 'int_n_bases_padding_around_interval' base pair away from the interval. 
        This functionality will allow analyzing all reads aligned to a single interval in the 'dict_it' to be analyzed in a single thread, thus abolishing the need to combine results from all threads to produce interval-level output, thus simplifying  downstream analysis. if chromosome name contains 'chr' prefix, please set 'flag_genome_without_chr_prefix_was_used_for_l_dict_it' to False
    'flag_genome_without_chr_prefix_was_used_for_l_dict_it' : (default = True) Please set this flag to False if the chromosome names of the given 'l_dict_it' annotations have a 'chr' prefix.

    returns
    bookmark, n_sam_records
    bookmark's columns : l_col = [ 'id_refname_start', 'int_refstart_start', 'int_index_chunk_start', 'id_refname_end', 'int_refstart_end', 'int_index_chunk_end' ]
    """
    if isinstance( l_dict_it, dict ) : # if a single dict_it is given, wrap the dict_it inside a list
        l_dict_it = [ l_dict_it ]

    ''' retrieve sequence length of reference sequence used in the input BAM file from the BAM header '''
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        dict_header = samfile.header.to_dict( )
    dict_seqname_to_len_seq = pd.DataFrame( dict_header[ 'SQ' ] ).set_index( 'SN' ).LN.to_dict( )

    """ define internal functions """
    def __chromosome_name_remove_chr__( str_name_chrom ) :
        ''' remove 'chr' prefix from the chromosome name '''
        if 'chr' == str_name_chrom[ : 3 ] :
            return str_name_chrom[ 3 : ]
        else :
            return str_name_chrom

    def __l_bookmark_to_df_bookmark__( l_bookmark ) :
        ''' convert l_bookmark to df_bookmark for a single chunk '''
        l_col = [ 'id_refname_start', 'int_refstart_start', 'int_index_chunk_start', 'id_refname_end', 'int_refstart_end', 'int_index_chunk_end' ]
        l_bookmark = np.array( l_bookmark, dtype = object ).reshape( int( len( l_bookmark ) / 2 ), 6 ) # reshape bookmark so that each record contains bookmark of the start point and end points
        df_bookmark = pd.DataFrame( l_bookmark, columns = l_col )
        return df_bookmark

    """ Compose list of bookmarks for iterating input BAM file using multiple threads """
    """ retrieve the positions and the id_reads at the start and stop points ('bookmark') for multiprocessing """
    ''' retrieve bookmarks to mark the start and ends of contigs for each chunnk '''
    l_bookmark = [ ] # empty bookmark
    flag_has_bookmark_been_initialized = False
    int_start_possible_boundary = None # start position of the possible boundary # will be set to non-None value when the desired number of records in the current chunk has been reached its target and the boundary is being set or being adjusted
    int_n_sam_record_count = 0
    int_index_chunk = 0
    id_refname_current = None # default reference name

    # iterate the entire BAM file to retrieve the total number of SAM records
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            ''' filter alignments based on mapq and flag '''
            if r.mapq < int_min_mapq_unique_mapped : # skip read whose mapq is below 'int_min_mapq_unique_mapped'
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
                    int_actual_n_sam_records_per_chunk = int_n_sam_record_count - 1
                    if verbose :
                        print( f'new chunk initiated (chunk {int_index_chunk + 1}) (chunk {int_index_chunk} has {int_n_sam_record_count} number of reads)' )
                    int_index_chunk += 1
                    int_n_sam_record_count = 1
                    bookmark = [ refname, 0, int_index_chunk ]
                    l_bookmark.append( bookmark ) 
                    ''' yield current bookmark and initialize the next bookmark '''
                    yield __l_bookmark_to_df_bookmark__( l_bookmark ), int_actual_n_sam_records_per_chunk
                    l_bookmark = [ ]
                    l_bookmark.append( bookmark )
                else :
                    bookmark = [ refname, 0, int_index_chunk ]
                    l_bookmark.append( bookmark ) 
                    l_bookmark.append( bookmark )

            """
            Search boundaries 
            """
            # initialize the next chunk if the 'int_n_sam_record_count' reached the assigned number of sam records for the current chunk (process)
            if int_n_sam_record_count > int_num_sam_records_for_each_chunk :
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
                refname_for_l_dict_it = __chromosome_name_remove_chr__( refname ) if flag_genome_without_chr_prefix_was_used_for_l_dict_it else refname # handle if BAM file's genome chromosome name convention and 'l_dict_it' annotations' chromosome name convention is different # retrieve refname without 'chr' prefix if 'flag_genome_without_chr_prefix_was_used_for_l_dict_it' is True

                for dict_it in l_dict_it :
                    if refname_for_l_dict_it in dict_it :
                        for start, end, anno_name in dict_it[ refname_for_l_dict_it ][ refstart ] :
                            l_adjusted_boundary_based_on_overlapped_interval.append( min( end + int_n_bases_padding_around_interval, int_len_seq_current_contig - 1 ) ) # maximum boundary is the end position of the contig
                            if verbose :
                                print( f"for chunk {int_index_chunk}, the possible boundary currently overlaps with {refname_for_l_dict_it}:{start}-{end}|{anno_name}" )
                # if there is any overlaps with the intervals, start adjusting boundaries
                if len( l_adjusted_boundary_based_on_overlapped_interval ) != 0 :
                    int_start_possible_boundary_previous = int_start_possible_boundary
                    # set the next possible boundary
                    int_start_possible_boundary = max( l_adjusted_boundary_based_on_overlapped_interval ) # use the farthest interval ends when setting the next potential boundary 
                    if verbose :
                        print( f'for chunk {int_index_chunk}, the possible boundary will change from {refname_for_l_dict_it}:{int_start_possible_boundary_previous} to {refname_for_l_dict_it}:{int_start_possible_boundary} ({int_start_possible_boundary - refstart}bp shift) when int_n_sam_record_count = {int_n_sam_record_count}' )
                    continue

                """
                Finalize boundary
                """
                ''' if the current position does not overlap with any intervals, set the current position as the boundary between chunks '''
                int_start_possible_boundary = None # set 'int_start_possible_boundary' to indicate boundary is not being searched
                int_actual_n_sam_records_per_chunk = int_n_sam_record_count - 1
                if verbose :
                    print( f'new chunk initiated (chunk {int_index_chunk + 1}) (chunk {int_index_chunk} has {int_n_sam_record_count} number of reads)' )
                # initialized the next chunk
                int_n_sam_record_count = 1
                int_index_chunk += 1
                # add a bookmark of the current position to mark the start of the new chunk
                bookmark = [ refname, refstart, int_index_chunk ]
                l_bookmark.append( bookmark )
                ''' yield current bookmark and initialize the next bookmark '''
                yield __l_bookmark_to_df_bookmark__( l_bookmark ), int_actual_n_sam_records_per_chunk
                l_bookmark = [ ]
                l_bookmark.append( bookmark )
    ''' yield final bookmark '''
    l_bookmark.append( [ '', -1, int_index_chunk ] )
    yield __l_bookmark_to_df_bookmark__( l_bookmark ), int_n_sam_record_count # yield the bookmark for the final chunk

def Call_Variant( 
    r, 
    dict_fasta_genome, 
    function_for_processing_reference_name = None
) :
    ''' # 2023-01-07 18:58:43 
    perform vcf-style variant calling of a single aligned read (pysam read object) using a given genome
    return a list of mutations and corrected read sequence
    
    'r' : AlignedSegment object returned by pysam
    'dict_fasta_genome' : dictionary of genome sequences
    'function_for_processing_reference_name' : a function that can be applied to reference name to make reference_name used in r consistent with those used in 'dict_fasta_genome' and the id_mut that will be returned by the current function. (e.g. a function that can remove 'chr' prefix from the reference name, if chromosome name without the 'chr' prefix is being used)

    --- returns ---
    'l_mut' : list of variants with the following nomenclature
    {refname}:{refpos}:{ref}>{alt}, where refpos, ref, alt follows nomenclature of VCF file
    '''
    pos_read, pos_ref = 0, 0 # current position in the extracted reference sequence and the current read
    
    # parse a sam record
    arr_qual = r.query_qualities
    seq = r.seq
    ref_name = r.reference_name
    if function_for_processing_reference_name is not None :
        ref_name = function_for_processing_reference_name( ref_name )
    ref_start = r.reference_start
    cigartuples = r.cigartuples
    alen = r.alen
    
    # retrieve a part of the reference sequence where the current read was aligned
    seq_ref = dict_fasta_genome[ ref_name ][ ref_start : ref_start + alen ] 
    
    l_mut = [ ]

    """
    # define interger representation of the CIGAR operations used in BAM files
    
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    """
    int_cigarop_M = 0
    int_cigarop_I = 1
    int_cigarop_D = 2
    int_cigarop_N = 3
    int_cigarop_S = 4
    int_cigarop_H = 5
    int_cigarop_P = 6
    int_cigarop_equal = 7
    int_cigarop_X = 8

    ns = dict( ) # create a namespace
    # initialilze the namespace
    ns[ 'pos_ref_variant_start' ] = None # 0-based coordinate of the start of the alternative allele on the reference
    ns[ 'pos_ref_variant_end' ] = None # 0-based coordinate of the end of the alternative allele on the reference
    ns[ 'alt' ] = ''
    def __update_alt( pos_ref : int, len_bases_ref : int = 0, bases_alt : str = '' ) :
        """ # 2023-01-07 17:51:14 
        """
        # initialize variant
        if ns[ 'pos_ref_variant_start' ] is None :
            ns[ 'pos_ref_variant_start' ] = pos_ref
            ns[ 'pos_ref_variant_end' ] = pos_ref
            ns[ 'alt' ] = ''
        ns[ 'pos_ref_variant_end' ] += len_bases_ref
        ns[ 'alt' ] += bases_alt
    def __flush_alt( ) :
        """ # 2023-01-07 17:24:20 
        flush variant allele
        """
        if ns[ 'pos_ref_variant_start' ] is not None :
            if len( ns[ 'alt' ] ) == 0 or ns[ 'pos_ref_variant_end' ] == ns[ 'pos_ref_variant_start' ] : # for simple insertion/deletion variants, add a single base pair of the reference before insertion or deletion to record the variant
                ns[ 'pos_ref_variant_start' ] -= 1
                ns[ 'alt' ] = seq_ref[ ns[ 'pos_ref_variant_start' ] ] + ns[ 'alt' ]
            l_mut.append( f"{ref_name}:{ref_start + 1 + ns[ 'pos_ref_variant_start' ]}:{seq_ref[ ns[ 'pos_ref_variant_start' ] : ns[ 'pos_ref_variant_end' ] ]}>{ns[ 'alt' ]}" ) # append the variant
            ns[ 'pos_ref_variant_start' ] = None
            ns[ 'pos_ref_variant_end' ] = None
            ns[ 'alt' ] = ''
    for int_oper, n_bases in cigartuples :
        if int_oper == int_cigarop_M :
            for _ in range( n_bases ) :
                str_base_read, str_base_ref, str_qual_read = seq[ pos_read ], seq_ref[ pos_ref ], arr_qual[ pos_read ]
                if str_base_read == str_base_ref :
                    __flush_alt( ) # flush a variant
                else :
                    __update_alt( pos_ref, len_bases_ref = 1, bases_alt = str_base_read )
                pos_ref += 1
                pos_read += 1
        elif int_oper == int_cigarop_N :
            __flush_alt( ) # flush a variant
            pos_ref += n_bases
        elif int_oper == int_cigarop_S :
            __flush_alt( ) # flush a variant
            pos_read += n_bases
        elif int_oper == int_cigarop_I :
            __update_alt( pos_ref, bases_alt = seq[ pos_read : pos_read + n_bases ] )
            pos_read += n_bases
        elif int_oper == int_cigarop_D :
            __update_alt( pos_ref, len_bases_ref = n_bases )
            pos_ref += n_bases
        elif int_oper == int_cigarop_H :
            __flush_alt( ) # flush a variant
            pass
        elif int_oper == int_cigarop_P :
            __flush_alt( ) # flush a variant
            pass
        elif int_oper == int_cigarop_equal :
            __flush_alt( ) # flush a variant
            pos_ref += n_bases
            pos_read += n_bases
        elif int_oper == int_cigarop_X :
            for _ in range( n_bases ) :
                str_base_read, str_qual_read = seq[ pos_read ], arr_qual[ pos_read ]
                __update_alt( pos_ref, 1, str_base_read )
                pos_ref += 1
                pos_read += 1
    __flush_alt( ) # flush a variant
    return l_mut
    
def Retrieve_Variant( r, dict_fasta_genome, set_mut_filtered = None, return_corrected_read_sequence = False, flag_ignore_indel = False, flag_return_as_string = True, flag_return_matched = False, function_for_processing_reference_name = None ) :
    ''' # 2022-01-15 00:37:28 
    perform variant calling of a single aligned read (pysam read object) using a given genome
    return a list of mutations and corrected read sequence
    
    'r' : AlignedSegment object returned by pysam
    'dict_fasta_genome' : dictionary of genome sequences
    'set_mut_filtered' : only consider mutations in the given 'set_mut_filtered'. only valid when 'return_corrected_read_sequence' = True
    'return_corrected_read_sequence' : return aligned read sequence after hard-clipping and filtering variants using 'set_mut_filtered'
    'flag_ignore_indel' : ignore insertion and deletion mutation calls
    'flag_return_as_string' : return a called mutation as a string in the following format: f"{reference_name}:{start}_{cigar_operation}_{mutation}". If set to False, return the called mutation as a tuple (reference_name, start, int_operation, mutation, quality)
    'flag_return_matched' : return 'id_mut' even when there is no mutation. the 'int_operation' for matched base is int_cigarop_M
    'function_for_processing_reference_name' : a function that can be applied to reference name to make reference_name used in r consistent with those used in 'dict_fasta_genome' and the id_mut that will be returned by the current function. (e.g. a function that can remove 'chr' prefix from the reference name)

    --- returns ---
    'l_mut' : list of variants with the following nomenclature
    # substitution:
    {refname}:{refpos}_{refbase}>{mutbase}

    # insertion:
    {refname}:{refpos_after_the_insertion}_ins_{inserted_bases}

    # deletion:
    {refname}:{refpos_of_the_first_base_of_the_deletion}_del_{number_of_deleted_bases}
    '''
    pos_read, pos_ref = 0, 0 # current position in the extracted reference sequence and the current read
    
    # parse a sam record
    arr_qual = r.query_qualities
    seq = r.seq
    ref_name = r.reference_name
    if function_for_processing_reference_name is not None :
        ref_name = function_for_processing_reference_name( ref_name )
    ref_start = r.reference_start
    cigartuples = r.cigartuples
    alen = r.alen
    
    # retrieve a part of the reference sequence where the current read was aligned
    seq_ref = dict_fasta_genome[ ref_name ][ ref_start : ref_start + alen ] 
    
    l_mut = [ ]

    """
    # define interger representation of the CIGAR operations used in BAM files
    
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    """
    int_cigarop_M = 0
    int_cigarop_I = 1
    int_cigarop_D = 2
    int_cigarop_N = 3
    int_cigarop_S = 4
    int_cigarop_H = 5
    int_cigarop_P = 6
    int_cigarop_equal = 7
    int_cigarop_X = 8

    if return_corrected_read_sequence : # correction mode
        str_seq_corrected_read = '' # corrected read sequence after hard clipping and filtering variants
        for int_oper, n_bases in cigartuples :
            if int_oper == int_cigarop_M :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = seq[ pos_read ], seq_ref[ pos_ref ], arr_qual[ pos_read ]
                    str_base_read_corrected = str_base_ref # default corrected read base = ref
                    if str_base_read != str_base_ref :
                        id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'X', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        if set_mut_filtered is None or id_mut in set_mut_filtered :
                            l_mut.append( id_mut )
                            str_base_read_corrected = str_base_read
                    elif flag_return_matched :
                        id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        l_mut.append( id_mut )
                    str_seq_corrected_read += str_base_read_corrected # add corrected read 
                    pos_ref += 1
                    pos_read += 1
            elif int_oper == int_cigarop_N :
                pos_ref += n_bases
            elif int_oper == int_cigarop_S :
                pos_read += n_bases
            elif int_oper == int_cigarop_I :
                id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_ins_{seq[ pos_read : pos_read + n_bases ]}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'I', seq[ pos_read : pos_read + n_bases ], arr_qual[ pos_read : pos_read + n_bases ] ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and ( set_mut_filtered is None or id_mut in set_mut_filtered ) :
                    l_mut.append( id_mut ) 
                    str_seq_corrected_read += seq[ pos_read : pos_read + n_bases ] # add inserted sequence (if insertion is valid)
                pos_read += n_bases
            elif int_oper == int_cigarop_D :
                id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_del_{n_bases}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'D', n_bases, None ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and ( set_mut_filtered is None or id_mut in set_mut_filtered ) :
                    l_mut.append( id_mut ) 
                else :
                    str_seq_corrected_read += seq_ref[ pos_ref : pos_ref + n_bases ] # add deleted reference sequence (if deletion is invalid)
                pos_ref += n_bases
            elif int_oper == int_cigarop_H :
                pass
            elif int_oper == int_cigarop_P :
                pass
            elif int_oper == int_cigarop_equal :
                if flag_return_matched :
                    str_base_read, str_qual_read = seq[ pos_read ], arr_qual[ pos_read ]
                    id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append( id_mut )
                str_seq_corrected_read += seq[ pos_read : pos_read + n_bases ] # add read sequences
                pos_ref += n_bases
                pos_read += n_bases
            elif int_oper == int_cigarop_X :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = seq[ pos_read ], seq_ref[ pos_ref ], arr_qual[ pos_read ]
                    str_base_read_corrected = str_base_ref # default corrected read base = ref
                    id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'X', str_base_read, str_qual_read ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    if set_mut_filtered is None or id_mut in set_mut_filtered :
                        l_mut.append( id_mut ) 
                        str_base_read_corrected = str_base_read
                    str_seq_corrected_read += str_base_read_corrected # add corrected read 
                    pos_ref += 1
                    pos_read += 1
        return l_mut, str_seq_corrected_read
    else : # normal mutation calling mode
        for int_oper, n_bases in cigartuples :
            if int_oper == int_cigarop_M :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = seq[ pos_read ], seq_ref[ pos_ref ], arr_qual[ pos_read ]
                    if str_base_read != str_base_ref :
                        l_mut.append( f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'X', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    elif flag_return_matched :
                        l_mut.append( f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'M', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
            elif int_oper == int_cigarop_N :
                pos_ref += n_bases
            elif int_oper == int_cigarop_S :
                pos_read += n_bases
            elif int_oper == int_cigarop_I :
                if not flag_ignore_indel :
                    l_mut.append( f"{ref_name}:{ref_start + 1 + pos_ref}_ins_{seq[ pos_read : pos_read + n_bases ]}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'I', seq[ pos_read : pos_read + n_bases ], arr_qual[ pos_read : pos_read + n_bases ] ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_read += n_bases
            elif int_oper == int_cigarop_D :
                if not flag_ignore_indel :
                    l_mut.append( f"{ref_name}:{ref_start + 1 + pos_ref}_del_{n_bases}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'D', n_bases, None ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_ref += n_bases
            elif int_oper == int_cigarop_H :
                pass
            elif int_oper == int_cigarop_P :
                pass
            elif int_oper == int_cigarop_equal :
                if flag_return_matched :
                    str_base_read, str_qual_read = seq[ pos_read ], arr_qual[ pos_read ]
                    id_mut = f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append( id_mut )
                pos_ref += n_bases
                pos_read += n_bases
            elif int_oper == int_cigarop_X :
                for _ in range( n_bases ) :
                    str_base_read, str_qual_read = seq[ pos_read ], arr_qual[ pos_read ]
                    l_mut.append( f"{ref_name}:{ref_start + 1 + pos_ref}_{str_base_ref}>{str_base_read}" if flag_return_as_string else ( ref_name, ref_start + pos_ref, 'X', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
        return l_mut

def Get_dict_pysam_alignmentheader_from_dict_seqname_to_len_seq( dict_seqname_to_len_seq ) :
    ''' # 2022-06-10 17:55:59 
    compose a dictionary that can be used to build an alignment header from 'dict_seqname_to_len_seq'
    '''
    dict_header = dict( )
    dict_header[ 'HD' ] = { 'VN': '1.6', 'SO': 'coordinate' }
    dict_header[ 'SQ' ] = list( { 'SN': seqname, 'LN': dict_seqname_to_len_seq[ seqname ] } for seqname in dict_seqname_to_len_seq )
    return dict_header

def Convert_mappy_alignment_to_pysam_aligned_segment( qname, seq, qual, mappy_alignment, samfile_header, start = None, flag_mark_as_supplementary_alignment = False ) :
    ''' # 2022-06-16 00:04:36 
    convert 'mappy_alignment' (mappy.Alignment record) to pysam.AlignedSegment object
    
    'flag_mark_as_supplementary_alignment' : mark the alignment record for supplementary alignment. for mappy, unlike minimap2 output, supplementary alignment is not marked by default (https://github.com/lh3/minimap2/issues/342). The first alignment can be considered primary, and all other alignment can be marked as supplementary.
    'start' : if a subset sequence of the 'query' sequence was searched using mappy, set this argument as the 0-based coordinate of the start of the subsequence
    '''
    len_seq = len( seq )
    # set default arguments
    if start is None :
        start = 0

    flag = 0 # initialize flag
    if not mappy_alignment.is_primary :
        flag ^= 1 << 8
    if flag_mark_as_supplementary_alignment :
        flag ^= 1 << 11
    if mappy_alignment.strand == -1 : # mark reverse complemented reads
        flag ^= 1 << 4
        flag_reverse_complemented = True
        # reverse complement reads
        seq = SEQ.Reverse_Complement( seq )
        qual = qual[ : : -1 ]
        # retrieve numbers of soft-clipped bases
        int_num_soft_clipped_bases_left, int_num_soft_clipped_bases_right = len_seq - mappy_alignment.q_en - start, mappy_alignment.q_st + start
    else :
        # retrieve numbers of soft-clipped bases
        int_num_soft_clipped_bases_left, int_num_soft_clipped_bases_right = mappy_alignment.q_st + start, len_seq - mappy_alignment.q_en - start
        
        
    # retrieve a modify cigar string
    cigar_str = ( str( int_num_soft_clipped_bases_left ) + 'S' if int_num_soft_clipped_bases_left > 0 else '' ) + mappy_alignment.cigar_str + ( str( int_num_soft_clipped_bases_right ) + 'S' if int_num_soft_clipped_bases_right else '' )

    l_sam = [ qname, flag, mappy_alignment.ctg, mappy_alignment.r_st + 1, mappy_alignment.mapq, cigar_str, '*', 0, 0, seq, qual ] # compose a basic sam record using list
    # add tag
    for name_attr, type_val in zip( [ 'NM', 'MD', 'cs' ], [ 'i', 'A', 'A' ] ) :
        val = getattr( mappy_alignment, name_attr )
        if val != '' : # if the current tag contains a valid value
            l_sam.append( ''.join( [ name_attr, ':', type_val, ':', str( val ) ] ) ) # add a SAM tag
    # add transcript strand tag
    if hasattr( mappy_alignment, 'trans_strand' ) :
        l_sam.append( 'ts:A:' + { -1 : '-', 0 : '.', 1 : '+' }[ mappy_alignment.trans_strand ] ) # add a SAM tag
    
    return pysam.AlignedSegment.fromstring( '\t'.join( map( str, l_sam ) ), samfile_header ) # return pysam AlignedSegment object

def Visualize_alignment( r, am = None, seq_ref_aligned = None, flag_display = True ) :
    ''' # 2022-06-12 16:01:22 
    'am' : minimap2 index used for alignment 
    'r' : pysam.AlignedSegment record
    'seq_ref_aligned' : if 'am' is not given, use this sequence as the aligned portion of the reference sequence
    'flag_display' : display the aligned record if True. return two strings, 'str_aligned_seq_ref' and 'str_aligned_seq_query' if False
    '''
    # retrieve aligned portion of reference and query sequences
    if am is not None :
        seq_ref = am.seq( r.reference_name, r.reference_start, r.reference_end )
    else :
        seq_ref = seq_ref_aligned
    seq_query = r.seq

    # initialize
    l_aligned_seq_ref = [ ]
    l_aligned_seq_query = [ ]

    pos_ref, pos_query = 0, 0

    l_mask_consumes_query = [ True, True, False, False, True, False, False, True, True ]
    l_mask_consumes_ref = [ True, False, True, True, False, False, False, True, True ]
    l_flag_record_alignment = [ True, True, True, False, False, False, False, True, True ]
    l_fill_value = [ '', '-', '-', f'__splicing__', '__soft_clipping__', '', '', '', '' ]
    for int_cigarop, int_num_bases in r.cigartuples :
        flag_consumes_query = l_mask_consumes_query[ int_cigarop ]
        flag_consumes_ref = l_mask_consumes_ref[ int_cigarop ]
        flag_record_alignment = l_flag_record_alignment[ int_cigarop ]
        str_fill_value = l_fill_value[ int_cigarop ]

        # annotate alignment
        if flag_record_alignment :
            if flag_consumes_ref :
                l_aligned_seq_ref.append( seq_ref[ pos_ref : pos_ref + int_num_bases ] )
            else :
                l_aligned_seq_ref.append( str_fill_value * int_num_bases )
            if flag_consumes_query :
                l_aligned_seq_query.append( seq_query[ pos_query : pos_query + int_num_bases ] )
            else :
                l_aligned_seq_query.append( str_fill_value * int_num_bases )
        elif len( str_fill_value ) > 0 :
            l_aligned_seq_ref.append( str_fill_value )
            l_aligned_seq_query.append( str_fill_value )

        # update coordinate
        if flag_consumes_query :
            pos_query += int_num_bases
        if flag_consumes_ref :
            pos_ref += int_num_bases
    
    # compose visual representation of the alignment
    str_aligned_seq_ref = ''.join( l_aligned_seq_ref )
    str_aligned_seq_query = ''.join( l_aligned_seq_query )
    
    if flag_display :
        from IPython.display import display, HTML
        # display a long text strings with overflow
        display(HTML("<div style='overflow: scroll; width: 100%; min-width: 1px; font-family:consolas; white-space: nowrap; display: block;'>" +
                    str_aligned_seq_ref + '<br>' + str_aligned_seq_query +
                     "</div>"))
        
    else :
        return str_aligned_seq_ref, str_aligned_seq_query

def Plot_alignment( l_r, figsize = ( 10, 1 ), len_head = 100, verbose = True, min_dx_abs = 1e-2, color = 'random' ) :
    """ # 2022-06-12 19:06:43 
    plot a given list of pysam.AlignedSegment records of a single read aligned against the same reference.
    
    'len_head' : the length of the arrow head
    'min_dx_abs' : the min absolute length of arrow stem
    """
    qname = l_r[ 0 ].qname
    fig, ax = plt.subplots( 1, 2, figsize = figsize )
    for i, r in enumerate( l_r ) :
        # retrieve a random color 
        c = Get_random_hex_color( ) if color == 'random' else color
        
        # draw an alignment on the reference
        y = i
        ax[ 0 ].arrow( x = r.query_alignment_start if not r.is_reverse > 0 else len( r.seq ) - r.query_alignment_end, y = y, dx = max( r.query_alignment_end - r.query_alignment_start, min_dx_abs ), dy = 0, head_width = 1, head_length = len_head, fc = c, ec = 'k', length_includes_head = True ) # direction of the arrow on the query is always + direction!
        
        # draw an alignment on the reference
        
        for ref_start, ref_end in Retrive_List_of_Mapped_Segments( r.cigartuples, r.reference_start )[ 0 ] :
            dx_ref = ref_end - ref_start if not r.is_reverse else ref_start - ref_end
            if dx_ref == 0 :
                dx_ref = min_dx_abs if r.is_reverse else - min_dx_abs
            ax[ 1 ].arrow( x = ref_start if not r.is_reverse else ref_end, y = y, dx = dx_ref, dy = 0, head_width = 1, head_length = len_head, fc = c, ec = 'k', length_includes_head = True )
        if verbose :
            print( r.flag, qname, r.mapq, r.query_alignment_start, r.query_alignment_end, r.reference_start, r.reference_end, r.is_reverse, Retrive_List_of_Mapped_Segments( r.cigar, r.reference_start, flag_is_cigartuples_from_mappy = False ) )
    ax[ 0 ].set_title( f'query aligned regions for\n{qname}' )
    ax[ 1 ].set_title( f'reference aligned regions for\n{qname}' )
    