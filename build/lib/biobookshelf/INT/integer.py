import numpy as np

def Split( int_total, n_parts ) :
    """ # 2021-07-28 20:06:05 
    split 'int_total' into 'n_parts' number of integers as equally as possible.  """
    int_remaining = int_total
    l = [ ]
    for _ in range( n_parts - 1 ) :
        int_split = round( int_total / n_parts )
        int_remaining -= int_split
        l.append( int_split )
    l.append( int_remaining )
    return l