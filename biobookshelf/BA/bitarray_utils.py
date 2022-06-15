from bitarray import bitarray
import numpy as np

def Find( ba, val = 1 ) :
    ''' # 2022-06-12 21:09:22 
    search the given value in the bitarray iteratively and return the start position of the occurrences
    '''
    len_ba = len( ba )
    l_int_pos_occurrence = [ ]
    int_pos_start = 0
    while int_pos_start < len_ba :
        int_pos_occurrence = ba.find( val, int_pos_start + 1 )
        if int_pos_occurrence < 0 :
            break
        else :
            l_int_pos_occurrence.append( int_pos_occurrence )
            int_pos_start = int_pos_occurrence
    return l_int_pos_occurrence

def Find_Segment( ba, background = 1, flag_use_bitwise_operation = True ) :
    ''' # 2022-06-12 19:36:38 
    find segment of a given bitarray for the given background 'background'. for example, when background = 1, find segment of 0
    
    'flag_use_bitwise_operation' : use bitwise operation. useful when the length of bitarray is very long
    '''
    len_ba = len( ba )
    if flag_use_bitwise_operation : 
        ''' the implementation using bitwise operation '''
        mask = ( ba ^ ba >> 1 ) # find boundary of value change
        l_int_pos = Find( mask, 1 )
        flag_start_is_a_segment = ba[ 0 ] != background
        if flag_start_is_a_segment :
            l_int_pos = [ 0 ] + l_int_pos
        if len( l_int_pos ) % 2 != 0 :
            l_int_pos += [ len_ba ]
        return np.array( l_int_pos ).reshape( ( int( len( l_int_pos ) / 2 ), 2 ) )
    else :
        ''' the implementation using simple iteration '''
        int_pos = 0 
        int_seg_start = None
        l_seg = [ ]
        while int_pos < len_ba :
            if int_seg_start is None and ba[ int_pos ] != background :
                int_seg_start = int_pos
            elif int_seg_start is not None and ba[ int_pos ] == background :
                l_seg.append( [ int_seg_start, int_pos ] )
                int_seg_start = None
            int_pos += 1
        if int_seg_start is not None :
            l_seg.append( [ int_seg_start, len_ba ] )
        return l_seg
