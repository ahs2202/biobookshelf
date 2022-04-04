from biobookshelf.main import *
import pandas as pd
import numpy as np


def Generate_Kmer( seq, window_size ) :
    """ 
    # 2021-02-20 15:14:13 
    generate a list of Kmers from the sequence with the given window_size  """
    return list( seq[ i : i + window_size ] for i in range( 0, len( seq ) - window_size + 1, 1 ) )


# In[ ]:


def Reverse_Complement( seq ) :
    ''' # 2021-02-04 11:47:19 
    Return reverse complement of 'seq' '''
    dict_dna_complement = { "A" : 'T', "T" : 'A', "C" : 'G', "G" : 'C', "N" : 'N', "-" : '-' }
    return ''.join( list( dict_dna_complement[ base ] for base in seq ) )[ : : -1 ]


# In[ ]:


# 2020-05-29 10:53:13 
dict_NGS__aa_to_codon = { 'I' : [ 'ATT', 'ATC', 'ATA' ], 'L' : [ 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG' ], 'V' : [ 'GTT', 'GTC', 'GTA', 'GTG' ], 'F' : [ 'TTT', 'TTC' ], 'M' : [ 'ATG' ], 'C' : [ 'TGT', 'TGC' ], 'A' : [ 'GCT', 'GCC', 'GCA', 'GCG' ], 'G' : [ 'GGT', 'GGC', 'GGA', 'GGG' ], 'P' : [ 'CCT', 'CCC', 'CCA', 'CCG' ], 'T' : [ 'ACT', 'ACC', 'ACA', 'ACG' ], 'S' : [ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ], 'Y' : [ 'TAT', 'TAC' ], 'W' : [ 'TGG' ], 'Q' : [ 'CAA', 'CAG' ], 'N' : [ 'AAT', 'AAC' ], 'H' : [ 'CAT', 'CAC' ], 'E' : [ 'GAA', 'GAG' ], 'D' : [ 'GAT', 'GAC' ], 'K' : [ 'AAA', 'AAG' ], 'R' : [ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ], '*' : [ 'TAA', 'TAG', 'TGA' ] }
dict_NGS__codon_to_aa = dict( )
for aa in dict_NGS__aa_to_codon :
    l_codon = dict_NGS__aa_to_codon[ aa ]
    for codon in l_codon : dict_NGS__codon_to_aa[ codon ] = aa
        
def Translate( seq, start = 0, append_stop_codon_to_returned_seq = False, print_message = True, return_error_value_if_stop_codon_not_found = True ) :
    ''' Translate nucleotide sequence until the given seq is end or stop codons are encountered.
    start = 0 : start position for potential protein sequence (for example, start of ATG start codon.) '''
    if type( seq ) is float : return np.nan # if input is invalid, return np.nan
    length_seq = len( seq ) # it seems the time complexity of the len(str) method is O(1)
    l_aa = list( )
    while True :
        if start + 3 > length_seq : 
            if print_message : print( 'Stop Codon Not Found' )
            if return_error_value_if_stop_codon_not_found : return np.nan
            else : break
        try : aa = dict_NGS__codon_to_aa[ seq[ start : start + 3 ] ] # try to translat the sequence
        except : 
            if return_error_value_if_stop_codon_not_found : return np.nan # if an error occur during translating and stop codon was not found, return np.nan if 'return_error_value_if_stop_codon_not_found'
            else : break
        if aa == '*' :
            if append_stop_codon_to_returned_seq : l_aa.append( aa )
            break
        else : l_aa.append( aa )
        start += 3
    return ''.join( l_aa )


# In[ ]:


# 2020-05-29 14:28:22 
def ORF_Find_All_Methionine_ORFs( seq, return_list = False ) :
    ''' Find all ORFs of a given sequence, and return it as a dataframe ( columns = [ 'start', 'end', 'translated_protein_seq' ], use 0-based coordinate for start and end ). 
    return found ORFs as 2D list if 'return_list' is set to True. '''
    l_l_value = [ [ 0, 0, np.nan ] ] # 2D array with a dummy starting ORF
    for int_atg_start_pos in STR.Find_all( seq, 'ATG' ) : # find all positions of subsequences starting with ATG (methionine start codon)
        prot_seq = Translate( seq, int_atg_start_pos, append_stop_codon_to_returned_seq = True, return_error_value_if_stop_codon_not_found = False, print_message = False )
        end_pos = int_atg_start_pos + 3 * len( prot_seq ) # retrive end position of CDS, including a stop codon if detected.
        if l_l_value[ - 1 ][ 1 ] == end_pos : continue # if current ORF shared the same stop position (and with downstream start position), do not include the ORF in the output, since current ORF is a fragment of the previous ORF
        else : l_l_value.append( [ int_atg_start_pos, end_pos, prot_seq ] ) # use 0-based coordinate
    return l_l_value[ 1 : ] if return_list else pd.DataFrame( l_l_value[ 1 : ], columns = [ 'start', 'end', 'translated_protein_seq' ] ) # return found ORFs without the dummy ORF at the start of the 2D array.
# 2020-05-29 14:28:22 
def ORF_Find_All_Methionine_ORFs_on_Both_Strands( seq, use_1_based_coordinate = True ) :
    ''' Find all ORFs start with methionine on both strand of the given sequence by using 'ORF_Find_All_Methionine_ORFs'. 
    use 1-based-coordinate for start and end positions if 'use_1_based_coordinate' is set to True '''
    df_plus = ORF_Find_All_Methionine_ORFs( seq ) # find all ORFs in the sequence and the reverse complement of the sequence
    df_minus = ORF_Find_All_Methionine_ORFs( Reverse_Complement( seq ) )
    df_plus[ 'strand' ] = '+'
    df_minus[ 'strand' ] = '-'
    int_length_seq = len( seq ) # convert coordinate of reverse complement of a sequence into that of the sequence
    s_start = int_length_seq - df_minus.end
    s_end = int_length_seq - df_minus.start
    df_minus[ 'start' ] = s_start
    df_minus[ 'end' ] = s_end

    df_orf = pd.concat( [ df_plus, df_minus ], ignore_index = True ) # combine ORFs from plus strand and the minus strand
    df_orf.sort_values( 'start', ignore_index = True, inplace = True )
    if use_1_based_coordinate : df_orf.start += 1 # use 1-based-coordinate for start and end positions if 'use_1_based_coordinate' is set to True
    return df_orf


# In[ ]:


def Trim_PolyA( seq, from_3_prime_end = True, return_length_of_polyA = False, str_repeated_base = 'A', int_lookup_window = 3 ) :
    ''' Trim PolyA sequence from the 3' end of a given sequence (DNA sequence in upper characters) if 'from_3_prime_end' is True or trim polyA from 5' end of the sequence if 'from_3_prime_end' is False
    if the program encounter base other than 'A', lookup next 'int_lookup_window' number of bases and check weather they are consecutive 'A' bases. '''
    if len( seq ) == 0 : # if the given seq is empty 
        if return_length_of_polyA : return 0
        else : return ''
    str_repeated_bases_in_lookup_window = str_repeated_base * int_lookup_window
    seq_polyA_at_5_prime_end = seq[ : : -1 ] if from_3_prime_end else seq
    for index, str_base in enumerate( seq_polyA_at_5_prime_end ) :
        if str_base != str_repeated_base :
            if seq_polyA_at_5_prime_end[ index + 1 : index + 1 + int_lookup_window ] == str_repeated_bases_in_lookup_window : continue # even though sequence is not 'A', look up the next bases, and if the next consequtive three bases are 'A', consider the character as heterogenity of PolyA and countinue counting number of 'A's in PolyA 
            else : break
    if return_length_of_polyA : return index
    else : return ( seq[  : - index ] if index > 0 else seq ) if from_3_prime_end else seq[ index : ]


# In[ ]:

def Detect_PolyT_Length( seq_after_softclipping, int_len_window_internal_polyT = 30, int_len_sliding_window_internal_polyT = 10, float_min_T_fraction = 0.8 ) :
    ''' # 2021-08-24 00:59:00 
    Detect polyT length of sequencing reads from an internal polyA priming event using a sliding window of a given length.
    '''
    ba = bitarray( len( seq_after_softclipping ) )
    ba.setall( 0 )

    for index, base in enumerate( seq_after_softclipping ) :
        ba[ index ] = base == 'T'

    int_len_internal_polyT = 0
    if ba[ : int_len_sliding_window_internal_polyT ].count( ) / int_len_sliding_window_internal_polyT >= float_min_T_fraction :
        int_len_internal_polyT = int_len_sliding_window_internal_polyT
        for index in range( 1, int_len_window_internal_polyT - int_len_sliding_window_internal_polyT + 1 ) :
            if ba[ index : index + int_len_sliding_window_internal_polyT ].count( ) / int_len_sliding_window_internal_polyT < float_min_T_fraction :
                break
            int_len_internal_polyT += 1
    return int_len_internal_polyT

dict_NGS_encoding_seq_to_int = { '' : 0,  'A' : 1 , 'C' : 2 , 'G' : 3 , 'T' : 4  }
def Encode_to_integer( seq, reverse_order_during_encoding = True ) :
    ''' convert sequence string ( should be in upper case ) with maximum length 27bp into integer (64 bit integer) in base-5 numeral system ('A' = 1, 'C' = 2, 'G' = 3, 'T' = 4). For algorithmical reason, encoded sequence is in an inverse order. '''
    int_encoded_seq = 0
    if not reverse_order_during_encoding : seq = seq[ : : -1 ] 
    for index, base in enumerate( seq ) :
        int_encoded_seq += dict_NGS_encoding_seq_to_int[ base ] * 5 ** index
    return int_encoded_seq

dict_NGS_decoding_int_to_seq = { 0 : '', 1 : 'A', 2 : 'C', 3 : 'G', 4 : 'T' }

def Decode_from_integer( int_encoded_seq, reverse_order_during_encoding = False ) :
    ''' convert integer (64 bit integer) in base-5 numeral system ('A' = 1, 'C' = 2, 'G' = 3, 'T' = 4) into sequence string ( in upper case ). For algorithmical reason, decoded sequence is in an inverse order. '''
    l_decoded_seq = list( )
    while True :
        decoded_base = dict_NGS_decoding_int_to_seq[ int_encoded_seq % 5 ]
        if not decoded_base : break
        l_decoded_seq.append( decoded_base )
        int_encoded_seq = int( int_encoded_seq / 5 )
    seq = ''.join( l_decoded_seq )
    return seq if reverse_order_during_encoding else seq[ : : -1 ]


# In[ ]:

## imported from biobookshelf.main

def COUNTER( l_values, dict_counter = None, ignore_float = True ) : # 2020-07-29 23:49:51 
    ''' Count values in l_values and return a dictionary containing count values. if 'dict_counter' is given, countinue counting by using the 'dict_counter'. if 'ignore_float' is True, ignore float values, including np.nan '''
    if dict_counter is None : dict_counter = dict( )
    if ignore_float : # if 'ignore_float' is True, ignore float values, including np.nan
        for value in l_values :
            if isinstance( value, float ) : continue # ignore float values
            if value in dict_counter : dict_counter[ value ] += 1
            else : dict_counter[ value ] = 1
    else : # faster counting by not checking type of value
        for value in l_values :
            if value in dict_counter : dict_counter[ value ] += 1
            else : dict_counter[ value ] = 1
    return dict_counter


def Calculate_Simple_Repeat_Proportion_in_a_read( seq, int_len_kmer = 4, int_kmer_count_threshold_for_repeat = 5 ) :
    ''' Calculate proportion of simple repeat in a read. simple repeat is detected by counting kmer of the given length and retrive kmer with count higher than the given threshold.  '''
    len_seq = len( seq )
    n_kmers = len_seq - int_len_kmer + 1
    dict_kmer_count = COUNTER( seq[ index : index + int_len_kmer ] for index in range( n_kmers ) ) # count kmers of the given length
    int_kmer_count_repeat = 0 # count kmers that have count values above the given threshold
    for kmer in dict_kmer_count :
        int_kmer_count = dict_kmer_count[ kmer ]
        if int_kmer_count > int_kmer_count_threshold_for_repeat : int_kmer_count_repeat += int_kmer_count
    float_ratio_of_repeat = int_kmer_count_repeat / n_kmers
    return float_ratio_of_repeat


# In[ ]:

def Cluster_with_Kmer( dict_seq_count, int_min_n_overlap_kmer, len_kmer, float_min_proportion_read_to_select_kmer_representative ) :
    """ # 2021-08-02 23:14:12 
    cluster given sequences (given as a dictionary containing counts of unique sequences) using kmer of a given length
    """
    dict_cluster = dict( )
    for seq in dict_seq_count :
        int_seq_count = dict_seq_count[ seq ]
        flag_assigned_to_cluster = False
        set_kmer = set( Generate_Kmer( seq, len_kmer ) )
        for name_cluster in dict_cluster :
            c = dict_cluster[ name_cluster ]
            n_overlap_kmer = len( c[ 'set_kmer' ].intersection( set_kmer ) )
            if int_min_n_overlap_kmer <= n_overlap_kmer :
                c[ 'seq_count' ][ seq ] = int_seq_count
                c[ 'n_seq' ] += int_seq_count
                n_seq = c[ 'n_seq' ]
                counter_kmer = COUNTER( list( set_kmer ) * int_seq_count, c[ 'counter_kmer' ] ) # update kmer count
                c[ 'counter_kmer' ] = counter_kmer
                c[ 'set_kmer' ] = set( kmer for kmer in counter_kmer if counter_kmer[ kmer ] / n_seq >= float_min_proportion_read_to_select_kmer_representative )
                flag_assigned_to_cluster = True
                break
        if not flag_assigned_to_cluster :
            c = dict( )
            c[ 'set_kmer' ] = set_kmer
            c[ 'n_seq' ] = int_seq_count
            c[ 'seq_count' ] = dict( )
            c[ 'seq_count' ][ seq ] = int_seq_count
            c[ 'counter_kmer' ] = COUNTER( list( set_kmer ) * int_seq_count ) 
            dict_cluster[ UUID( ) ] = c
    return dict_cluster

def Iterate_Kmer( seq, window_size, flag_return_start_and_end_positions = False, flag_generate_kmer_for_reverse_complement_too = False ) :
    """ 
    # 2021-09-18 11:00:58 
    return an interator generating Kmers from the sequence with the given window_size  
    
    returns
    ------------
    when 'flag_return_start_and_end_positions' = True, 
    returns kmer, pos_start, pos_end, flag_reverse_completed
    
    when 'flag_return_start_and_end_positions' = False,
    returns kmer, flag_reverse_completed
    """
    for i in range( 0, len( seq ) - window_size + 1, 1 ) :
        kmer = seq[ i : i + window_size ] 
        yield ( kmer, i, i + window_size, False ) if flag_return_start_and_end_positions else ( kmer, False )
    ''' if 'flag_generate_kmer_for_reverse_complement_too' is True, iterate reverse complement of the given sequence, too '''
    if flag_generate_kmer_for_reverse_complement_too :
        seq_rc = Reverse_Complement( seq )
        len_seq = len( seq )
        for i in range( 0, len_seq - window_size + 1, 1 ) :
            kmer = seq_rc[ i : i + window_size ] 
            yield ( kmer, len_seq - window_size - i, len_seq - i, True ) if flag_return_start_and_end_positions else ( kmer, True )