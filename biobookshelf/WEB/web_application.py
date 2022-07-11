from biobookshelf.main import *

def Round_Float( df, l_col_scientific_notations, l_col_typical_notation, n_significant_digits_scientific_notation = 3, n_significant_digits_typical_notation = 3, inplace = False ) :
    ''' # 2021-07-19 17:38:40 
    round float with a given number of significant digits and convert floating point numbers to strings '''
    if not inplace :
        df = deepcopy( df )
    str_format_scientific_notation = "{:." + str( int( n_significant_digits_scientific_notation ) ) + "e}"
    str_format_typical_notation = '{:.' + str( int( n_significant_digits_typical_notation ) ) + 'f}'
    for col in l_col_scientific_notations :
        df[ col ] = list( '' if np.isnan( value ) else str_format_scientific_notation.format( value ) for value in df[ col ].values )
    for col in l_col_typical_notation :
        df[ col ] = list( '' if np.isnan( value ) else str_format_typical_notation.format( value ) for value in df[ col ].values )
    return df

def Index_and_Base64_Encode( df_to_be_indexed, l_col_index, path_prefix_output, path_folder_temp = '/tmp/', flag_delete_temp_folder = True ) :
    """ # 2022-01-24 13:08:12 
    'df_to_be_indexed' : dataframe to be exported to base64 encoded file and indexed with values in 'l_col_index'
    'l_col_index' : list of columns for indexing 'df_to_be_indexed'
    'path_prefix_output' : directory prefix for an indexed base64 encoded file and a base64 encoded index file
    'path_folder_temp' : directory where a temporary folder will be created and removed
    'flag_delete_temp_folder' : flag indicating whether the temporary folder should be removed 
    """

    # retrieve path_folder_wd
    path_folder_temp = os.path.abspath( path_folder_temp )
    if path_folder_temp[ -1 ] != '/' :
        path_folder_temp += '/'
    str_uuid = UUID( )
    path_folder_temp = f"{path_folder_temp}{str_uuid}/"
    os.makedirs( path_folder_temp, exist_ok = True )

    df_to_be_indexed = deepcopy( df_to_be_indexed )
    # retrieve unique indices from values in the 'l_col_index'
    l_index = LIST_COUNT( df_to_be_indexed[ l_col_index ].values, duplicate_filter = None, sort_series_by_values = True ).index.values # sort by number of counts so that the most frequent access is close to the start of the file

    df_to_be_indexed.set_index( l_col_index, inplace = True )
    df_to_be_indexed.sort_index( inplace = True ) # sort by index for faster indexing
    
    # base64 encode each individual file
    for int_index, t_index in enumerate( l_index ) :
        df = df_to_be_indexed.loc[ t_index ]
        df.reset_index( drop = False, inplace = True )
        path_prefix_file = f"{path_folder_temp}{int_index}"
        df.T.to_csv( f"{path_prefix_file}.tsv.gz", sep = '\t', index = True, header = None ) # save transposed array (each 'row' is column, and the first element in each 'row' is the column name)
        Base64_Encode( f"{path_prefix_file}.tsv.gz", f"{path_prefix_file}.tsv.gz.base64.txt", header = ' ' )

    # retrieve file size of base64 encoded chucks
    df_file_base64 = GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_temp}*.tsv.gz.base64.txt", retrieve_file_size = True )
    df_file_base64.wildcard_0 = df_file_base64.wildcard_0.astype( int ) # retrieve integer index
    df_file_base64.sort_values( 'wildcard_0', inplace = True ) # sort by gene_name

    # concatanate base64 encoded files in the specified order
    OS_FILE_Combine_Files_in_order( df_file_base64.path.values, f"{path_prefix_output}.tsv.gz.base64.concatanated.txt", overwrite_existing_file = True )

    # write an index file describing the byte positions of each gene_symbol in the concatanated file
    int_byte_accumulated = 0
    l_l = [ ]
    flag_multiindex = len( l_col_index ) > 1 # retrieve flag indicating whether a multi-index was used.
    for int_index, size_in_bytes in df_file_base64[ [ 'wildcard_0', 'size_in_bytes' ] ].values :
        l_l.append( list( l_index[ int_index ] ) + [ int_byte_accumulated, int_byte_accumulated + size_in_bytes ] ) # index_byte uses 0-based coordinates
        int_byte_accumulated += size_in_bytes # update accumulated number of bytes
    df_index_byte = pd.DataFrame( l_l, columns = l_col_index + [ 'index_byte_start', 'index_byte_end' ] )

    df_index_byte.T.to_csv( f"{path_folder_temp}index.tsv.gz", sep = '\t', index = True, header = False )
    Base64_Encode( f"{path_folder_temp}index.tsv.gz", f"{path_prefix_output}.index.tsv.gz.base64.txt" ) # convert binary file into text using base64 encoding

    if flag_delete_temp_folder :
        # remove temporary folder
        shutil.rmtree( path_folder_temp )
        
def Index_Chunks_and_Base64_Encode( df_to_be_chunked_and_indexed, int_num_rows_for_each_chunk, path_prefix_output, flag_include_header = False, path_folder_temp = '/tmp/', flag_delete_temp_folder = True ) :
    """ # 2022-06-20 10:53:13 
    'df_to_be_chunked_and_indexed' : dataframe to be exported to base64 encoded file and indexed with values in 'l_col_index'
    'l_col_index' : list of columns for indexing 'df_to_be_chunked_and_indexed'
    'path_prefix_output' : directory prefix for an indexed base64 encoded file and a base64 encoded index file
    'flag_include_header' : include header in the output file
    'path_folder_temp' : directory where a temporary folder will be created and removed
    'flag_delete_temp_folder' : flag indicating whether the temporary folder should be removed 
    """

    # retrieve path_folder_wd
    path_folder_temp = os.path.abspath( path_folder_temp )
    if path_folder_temp[ -1 ] != '/' :
        path_folder_temp += '/'
    str_uuid = UUID( )
    path_folder_temp = f"{path_folder_temp}{str_uuid}/"
    os.makedirs( path_folder_temp, exist_ok = True )

    df_to_be_chunked_and_indexed = deepcopy( df_to_be_chunked_and_indexed )
    int_num_rows = len( df_to_be_chunked_and_indexed )

    # calculate the number of chunk
    int_num_chunk = int( np.ceil( int_num_rows / int_num_rows_for_each_chunk ) )
    # base64 encode each individual file
    for int_index in range( int_num_chunk ) :
        df = df_to_be_chunked_and_indexed.iloc[ int_index * int_num_rows_for_each_chunk : ( int_index + 1 ) * int_num_rows_for_each_chunk ]
        path_prefix_file = f"{path_folder_temp}{int_index}"
        df.T.to_csv( f"{path_prefix_file}.tsv.gz", sep = '\t', index = flag_include_header, header = None ) # save transposed array (each 'row' is column, and the first element in each 'row' is the column name) # to aid parsing in browser side
        Base64_Encode( f"{path_prefix_file}.tsv.gz", f"{path_prefix_file}.tsv.gz.base64.txt", header = ' ' ) # the additional ' ' space is to help separate file contents when using range queries

    # retrieve file size of base64 encoded chucks
    df_file_base64 = GLOB_Retrive_Strings_in_Wildcards( f"{path_folder_temp}*.tsv.gz.base64.txt", retrieve_file_size = True )
    df_file_base64.wildcard_0 = df_file_base64.wildcard_0.astype( int ) # retrieve integer index
    df_file_base64.sort_values( 'wildcard_0', inplace = True ) # sort by gene_name

    # concatanate base64 encoded files in the specified order
    OS_Run( [ 'cat' ] + list( df_file_base64.path.values ), path_file_stdout = f"{path_prefix_output}.tsv.gz.base64.concatenated.txt" )

    # write an index file describing the byte positions of each gene_symbol in the concatanated file
    int_byte_accumulated = 0
    l_l = [ ]
    for int_index, size_in_bytes in df_file_base64[ [ 'wildcard_0', 'size_in_bytes' ] ].values :
        l_l.append( [ int_index, int_byte_accumulated, int_byte_accumulated + size_in_bytes ] ) # index each chunk # index_byte uses 0-based coordinates
        int_byte_accumulated += size_in_bytes # update accumulated number of bytes
    df_index_byte = pd.DataFrame( l_l, columns = [ 'index_chunk', 'index_byte_start', 'index_byte_end' ] )

    df_index_byte.T.to_csv( f"{path_folder_temp}index.tsv.gz", sep = '\t', index = True, header = False )
    Base64_Encode( f"{path_folder_temp}index.tsv.gz", f"{path_prefix_output}.index.tsv.gz.base64.txt" ) # convert binary file into text using base64 encoding

    if flag_delete_temp_folder :
        # remove temporary folder
        shutil.rmtree( path_folder_temp )