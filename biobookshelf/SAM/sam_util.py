from biobookshelf.main import *

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
                
                
def Retrive_List_of_Mapped_Segments( cigartuples, pos_start, return_1_based_coordinate = False, flag_pos_start_0_based_coord = True, flag_return_splice_junction = False ) :
    ''' # 2021-09-07 10:17:55 
    return l_seq and int_total_aligned_length for given cigartuples (returned by pysam cigartuples) and 'pos_start' (0-based coordinates, assuming pos_start is 0-based coordinate)
    'return_1_based_coordinate' : return 1-based coordinate, assuming 'pos_start' is 0-based coordinate (pysam returns 0-based coordinate)
    'flag_return_splice_junction' : additionally return list of splice junction tuples
    '''
    if return_1_based_coordinate and flag_pos_start_0_based_coord : # 0-based > 1-based
        start -= 1
    l_seg, start, int_aligned_length, int_total_aligned_length = list( ), pos_start, 0, 0
    for operation, length in cigartuples :
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