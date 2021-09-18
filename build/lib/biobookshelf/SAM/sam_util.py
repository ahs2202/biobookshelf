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
                