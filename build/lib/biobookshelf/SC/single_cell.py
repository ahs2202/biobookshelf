from biobookshelf.main import *
import pandas as pd
import numpy as np

def MTX_10X_Read( dir_folder_mtx_10x, verbose = False ) :
    ''' # 2021-11-24 13:00:13 
    read 10x count matrix
    'dir_folder_mtx_10x' : a folder containing files for 10x count matrix
    return df_mtx, df_feature
    '''
    # handle inputs
    if dir_folder_mtx_10x[ -1 ] != '/' :
        dir_folder_mtx_10x += '/'
    
    # define input file directories
    dir_file_bc = f'{dir_folder_mtx_10x}barcodes.tsv.gz'
    dir_file_feature = f'{dir_folder_mtx_10x}features.tsv.gz'
    dir_file_mtx = f'{dir_folder_mtx_10x}matrix.mtx.gz'

    # check whether all required files are present
    if sum( list( not os.path.exists( dir_folder ) for dir_folder in [ dir_file_bc, dir_file_feature, dir_file_mtx ] ) ) :
        if verbose :
            print( f'required file(s) is not present at {dir_folder_mtx_10x}' )

    # read mtx file as a tabular format
    df_mtx = pd.read_csv( dir_file_mtx, sep = ' ', comment = '%' )
    df_mtx.columns = [ 'id_row', 'id_column', 'read_count' ]

    # read barcode and feature information
    df_bc = pd.read_csv( dir_file_bc, sep = '\t', header = None )
    df_bc.columns = [ 'barcode' ]
    df_feature = pd.read_csv( dir_file_feature, sep = '\t', header = None )
    df_feature.columns = [ 'id_feature', 'feature', '10X_type' ]

    # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx[ 'barcode' ] = df_mtx.id_column.apply( MAP.Map( DICTIONARY_Build_from_arr( df_bc.barcode.values, index_start = 1 ) ).a2b ) # mapping using 1 based coordinates (0->1 based coordinate )
    df_mtx[ 'id_feature' ] = df_mtx.id_row.apply( MAP.Map( DICTIONARY_Build_from_arr( df_feature.id_feature.values, index_start = 1 ) ).a2b ) 
    df_mtx.drop( columns = [ 'id_row', 'id_column' ], inplace = True ) # drop unnecessary columns
    
    return df_mtx, df_feature
def MTX_10X_Write( df_mtx, df_feature, dir_folder_output_mtx_10x ) :
    """ # 2021-11-24 12:57:30 
    'df_feature' should contains the following column names : [ 'id_feature', 'feature', '10X_type' ]
    'df_mtx' should contains the following column names : [ 'id_feature', 'barcode', 'read_count' ]
    'dir_folder_output_mtx_10x' : an output folder directory where the mtx_10x files will be written

    """
    df_mtx = deepcopy( df_mtx ) # create a copy of df_mtx before modification

    # create an output folder
    os.makedirs( dir_folder_output_mtx_10x, exist_ok = True )

    ''' save barcode file '''
    # retrieve list of barcodes
    arr_barcode = LIST_COUNT( df_mtx.barcode, duplicate_filter = None ).index.values
    pd.DataFrame( arr_barcode ).to_csv( f"{dir_folder_output_mtx_10x}barcodes.tsv.gz", sep = '\t', index = False, header = False ) 

    ''' save feature file '''
    # compose a feature dataframe
    df_feature[ [ 'id_feature', 'feature', '10X_type' ] ].to_csv( f"{dir_folder_output_mtx_10x}features.tsv.gz", sep = '\t', index = False, header = False ) # save as a file
    # retrieve list of features
    arr_id_feature = df_feature.id_feature.values

    ''' save matrix file '''
    # convert feature and barcode to integer indices
    df_mtx.id_feature = df_mtx.id_feature.apply( MAP.Map( DICTIONARY_Build_from_arr( arr_id_feature, order_index_entry = False ) ).a2b ) # 0-based coordinates
    df_mtx.barcode = df_mtx.barcode.apply( MAP.Map( DICTIONARY_Build_from_arr( arr_barcode, order_index_entry = False ) ).a2b ) # 0-based coordinates
    # save count matrix as a gzipped matrix market format
    row, col, data = df_mtx[ [ 'id_feature', 'barcode', 'read_count' ] ].values.T
    sm = scipy.sparse.coo_matrix( ( data, ( row, col ) ), shape = ( len( arr_id_feature ), len( arr_barcode ) ) )
    scipy.io.mmwrite( f"{dir_folder_output_mtx_10x}matrix", sm )
    # remove previous output file to overwrite the file
    dir_file_mtx_output = f"{dir_folder_output_mtx_10x}matrix.mtx.gz"
    if os.path.exists( dir_file_mtx_output ) :
        os.remove( dir_file_mtx_output )
    OS_Run( [ 'gzip', f"{dir_folder_output_mtx_10x}matrix.mtx" ] ) # gzip the mtx file
def SCANPY_Retrieve_Markers_as_DataFrame( adata ) :
    ''' # 2022-02-15 14:40:02 
    receive scanpy anndata and return a dataframe contianing marker genes 
    
    --- return --- 
    df_marker : a dataframe contianing marker genes 
    '''
    l_df = [ ]
    for index_clus, name_clus in enumerate( adata.uns["rank_genes_groups"]['names'].dtype.names ) : 
        df = pd.DataFrame( dict( ( name_col, adata.uns["rank_genes_groups"][ name_col ][ name_clus ] ) for name_col in ['logfoldchanges', 'names', 'pvals', 'pvals_adj', 'scores' ] ) )
        df[ "name_clus" ] = name_clus 
        df[ "index_clus" ] = index_clus
        l_df.append( df )
    df_marker = pd.concat( l_df )
    return df_marker
def __function_for_adjusting_thresholds_for_filtering_empty_droplets__( dir_folder_mtx_10x_output, min_counts, min_features, min_cells ) :
    ''' # 2022-02-23 14:26:07 
    This function is intended for the use in 'MTX_10X_Filter' function for filtering cells from the 10X dataset (before chromium X, 10,000 cells per channel)
    
    Assuming a typical number of droplets in a experiment is 100,000, adjust 'min_counts' to reduce the number of filtered cells below 'int_max_num_cells' 
    '''
    s_count = pd.read_csv( f"{dir_folder_mtx_10x_output}dict_id_column_to_count.before_filtering.tsv.gz", sep = '\t', header = None, index_col = 0 )[ 1 ].sort_values( ascending = False ).iloc[ : 100000 ]
    
    int_max_num_cells = 20000 # maximum number of allowed cells
    min_counts_maximum = 2000
    def function_for_increasing_min_counts( min_counts ) :
        return min_counts * 2
    while True :
        ''' increase threshold if the number of filtered cells is larger than 'int_max_num_cells' '''
        if len( s_count[ s_count > min_counts ] ) > int_max_num_cells and min_counts < min_counts_maximum :
            min_counts = function_for_increasing_min_counts( min_counts )
        else :
            break
    return min_counts, min_features, min_cells
def MTX_10X_Split( dir_folder_mtx_10x_output, int_max_num_entries_for_chunk = 10000000 ) :
    ''' # 2022-02-22 00:41:20 
    split input mtx file into multiple files and write a flag file indicating the splitting has been completed. 
    return the list of split mtx files
    '''
    dir_file_flag = f"{dir_folder_mtx_10x_output}matrix.mtx.gz.split.flag"
    if not os.path.exists( dir_file_flag ) : # check whether the flag exists
        index_mtx_10x = 0
        newfile = gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", 'wb' )
        l_dir_file_mtx_10x = [ f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz" ]
        int_num_entries_written_for_the_current_chunk = 0
        with gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz", 'rb' ) as file :
            while True :
                line = file.readline( ) # binary string
                if len( line ) == 0 :
                    newfile.close( ) # close the output file
                    break
                ''' write the line to the current chunk and update the number of entries written for the current chunk '''
                newfile.write( line )
                int_num_entries_written_for_the_current_chunk += 1
                ''' initialize the next chunk if a sufficient number of entries were written '''
                if int_num_entries_written_for_the_current_chunk >= int_max_num_entries_for_chunk :
                    newfile.close( ) # close the output file
                    index_mtx_10x += 1
                    int_num_entries_written_for_the_current_chunk = 0
                    newfile = gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", 'wb' )
                    l_dir_file_mtx_10x.append( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz" )
        with open( dir_file_flag, 'w' ) as file :
            file.write( 'completed' )
    else :
        ''' retrieve the list of split mtx files '''
        df = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.*.gz" )
        df.wildcard_0 = df.wildcard_0.astype( int )
        df.sort_values( 'wildcard_0', ascending = True, inplace = True )
        l_dir_file_mtx_10x = df.dir.values
    return l_dir_file_mtx_10x
dict_id_feature_to_index_feature = dict( )
def __MTX_10X_Combine__renumber_feature_mtx_10x__( dir_file_input, dir_folder_mtx_10x_output ) :
    '''
    internal function for MTX_10X_Combine
    # 2022-02-22 00:38:33 
    '''
#     dict_id_feature_to_index_feature = PICKLE_Read( f'{dir_folder_mtx_10x_output}dict_id_feature_to_index_feature.pickle' ) # retrieve id_feature to index_feature mapping 
    for dir_folder_mtx_10x, int_total_n_barcodes_of_previously_written_matrices, index_mtx_10x in pd.read_csv( dir_file_input, sep = '\t' ).values :
        # directly write matrix.mtx.gz file without header
        with gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", 'wb' ) as newfile :
            arr_id_feature = pd.read_csv( f'{dir_folder_mtx_10x}features.tsv.gz', sep = '\t', header = None ).values[ :, 0 ] # retrieve a list of id_feature for the current dataset
            with gzip.open( f'{dir_folder_mtx_10x}matrix.mtx.gz', 'rb' ) as file : # retrieve a list of features
                line = file.readline( ).decode( ) # read the first line
                # if the first line of the file contains a comment line, read all comment lines and a description line following the comments.
                if line[ 0 ] == '%' :
                    # read comment and the description line
                    while True :
                        if line[ 0 ] != '%' :
                            break
                        line = file.readline( ).decode( ) # read the next line
                    line = file.readline( ).decode( ) # discard the description line and read the next line
                # process entries
                while True :
                    if len( line ) == 0 :
                        break
                    index_row, index_col, int_value = tuple( map( int, line.strip( ).split( ) ) ) # parse each entry of the current matrix 
                    newfile.write( ( ' '.join( tuple( map( str, [ dict_id_feature_to_index_feature[ arr_id_feature[ index_row - 1 ] ], index_col + int_total_n_barcodes_of_previously_written_matrices, int_value ] ) ) ) + '\n' ).encode( ) ) # translate indices of the current matrix to that of the combined matrix            
                    line = file.readline( ).decode( ) # read the next line
def MTX_10X_Combine( dir_folder_mtx_10x_output, * l_dir_folder_mtx_10x_input, int_num_threads = 15, flag_split_mtx = True, int_max_num_entries_for_chunk = 10000000 ) :
    '''
    # 2022-02-22 00:38:36 
    Combine 10X count matrix files from the given list of folders and write combined output files to the given output folder 'dir_folder_mtx_10x_output'
    If there are no shared cells between matrix files, a low-memory mode will be used. The output files will be simply combined since no count summing operation is needed. Only feature matrix will be loaded and updated in the memory.
    'id_feature' should be unique across all features
    
    'int_num_threads' : number of threads to use when combining datasets. multiple threads will be utilized only when datasets does not share cells and thus can be safely concatanated.
    'flag_split_mtx' : split the resulting mtx file so that the contents in the output mtx file can be processed in parallel without ungzipping the mtx.gz file and spliting the file.
    '''
    
    # create an output folder
    os.makedirs( dir_folder_mtx_10x_output, exist_ok = True ) 

    """ retrieve cell barcodes of all 10X matrices and check whether cell barcodes are not shared between matrices """
    int_total_n_barcodes_of_previously_written_matrices = 0 # follow the number of barcodes that are previously written
    l_int_total_n_barcodes_of_previously_written_matrices = [ ] # calculate the number of barcodes of the previous dataset in the combined mtx.
    set_barcode = set( ) # update a set of unique barcodes
    for dir_folder_mtx_10x in l_dir_folder_mtx_10x_input :
        arr_barcode = pd.read_csv( f'{dir_folder_mtx_10x}barcodes.tsv.gz', sep = '\t', header = None ).squeeze( "columns" ).values # retrieve a list of features
        set_barcode.update( arr_barcode ) # update a set of barcodes
        l_int_total_n_barcodes_of_previously_written_matrices.append( int_total_n_barcodes_of_previously_written_matrices )
        int_total_n_barcodes_of_previously_written_matrices += len( arr_barcode ) # update the number of barcodes 

    ''' check whether there are shared cell barcodes between matrices and set a flag for entering a low-memory mode '''
    flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs = len( set_barcode ) == int_total_n_barcodes_of_previously_written_matrices # update flag

    if flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs :
        """ low-memory mode """
        """ write a combined barcodes.tsv.gz """
        OS_FILE_Combine_Files_in_order( list( f"{dir_folder_mtx_10x}barcodes.tsv.gz" for dir_folder_mtx_10x in l_dir_folder_mtx_10x_input ), f"{dir_folder_mtx_10x_output}barcodes.tsv.gz", overwrite_existing_file = True )

        ''' collect a set of unique features and a list of id_feature for each 10X matrix '''
        # assumes id_feature is unique for all features across 10X matrices
        set_t_feature = set( ) # update a set of feature tuples
        for dir_folder_mtx_10x in l_dir_folder_mtx_10x_input :
            arr_feature = pd.read_csv( f'{dir_folder_mtx_10x}features.tsv.gz', sep = '\t', header = None ).values # retrieve a list of features
            set_t_feature.update( list( map( tuple, arr_feature ) ) ) # update a set of feature tuples

        """ write a combined features.tsv.gz """
        l_t_feature = list( set_t_feature ) # convert set to list
        with gzip.open( f"{dir_folder_mtx_10x_output}features.tsv.gz", 'wb' ) as newfile :
            for t_feature in l_t_feature :
                newfile.write( ( '\t'.join( t_feature ) + '\n' ).encode( ) )

        """ build a mapping of id_feature to index_feature, which will be consistent across datasets """
        global dict_id_feature_to_index_feature # use global variable for multiprocessing
        dict_id_feature_to_index_feature = dict( ( t_feature[ 0 ], index_feature + 1 ) for index_feature, t_feature in enumerate( l_t_feature ) )
        PICKLE_Write( f'{dir_folder_mtx_10x_output}dict_id_feature_to_index_feature.pickle', dict_id_feature_to_index_feature ) # save id_feature to index_feature mapping as a pickle file

        ''' collect the number of entries for each 10X matrix '''
        int_total_n_entries = 0 
        for dir_folder_mtx_10x in l_dir_folder_mtx_10x_input :
            with gzip.open( f'{dir_folder_mtx_10x}matrix.mtx.gz', 'rb' ) as file : # retrieve a list of features
                file.readline( ), file.readline( )
                int_total_n_entries += int( file.readline( ).decode( ).strip( ).split( )[ 2 ] ) # update the total number of entries

        """ write a part of a combined matrix.mtx.gz for each dataset using multiple processes """
        # compose inputs for multiprocessing
        df_input = pd.DataFrame( { 'dir_folder_input_mtx_10x' : l_dir_folder_mtx_10x_input, 'int_total_n_barcodes_of_previously_written_matrices' : l_int_total_n_barcodes_of_previously_written_matrices, 'index_mtx_10x' : np.arange( len( l_int_total_n_barcodes_of_previously_written_matrices ) ) } )
        Multiprocessing( df_input, __MTX_10X_Combine__renumber_feature_mtx_10x__, int_num_threads, global_arguments = [ dir_folder_mtx_10x_output ] )
#         os.remove( f'{dir_folder_mtx_10x_output}dict_id_feature_to_index_feature.pickle' ) # remove pickle file
        
        """ combine parts and add the MTX file header to compose a combined mtx file """
        df_file = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.*.gz" )
        df_file.wildcard_0 = df_file.wildcard_0.astype( int )
        df_file.sort_values( 'wildcard_0', inplace = True )
        # if 'flag_split_mtx' is True, does not delete the split mtx files
        OS_FILE_Combine_Files_in_order( df_file.dir.values, f"{dir_folder_mtx_10x_output}matrix.mtx.gz", delete_input_files = not flag_split_mtx, header = f"%%MatrixMarket matrix coordinate integer general\n%\n{len( l_t_feature )} {len( set_barcode )} {int_total_n_entries}\n" ) # combine the output mtx files in the order
        # write a flag indicating that the current output directory contains split mtx files
        with open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.split.flag", 'w' ) as file :
            file.write( 'completed' )
        
    else :
        ''' normal operation mode perfoming count merging operations '''
        l_df_mtx, l_df_feature = [ ], [ ]
        for dir_folder_mtx_10x in l_dir_folder_mtx_10x_input :
            df_mtx, df_feature = MTX_10X_Read( dir_folder_mtx_10x )
            l_df_mtx.append( df_mtx ), l_df_feature.append( df_feature )

        # combine mtx
        df_mtx = pd.concat( l_df_mtx )
        df_mtx = df_mtx.groupby( [ 'barcode', 'id_feature' ] ).sum( )
        df_mtx.reset_index( drop = False, inplace = True )

        # combine features
        df_feature = pd.concat( l_df_feature )
        df_feature.drop_duplicates( inplace = True )

        MTX_10X_Write( df_mtx, df_feature, dir_folder_mtx_10x_output )
        
        # split a matrix file into multiple files
        MTX_10X_Split( dir_folder_mtx_10x_output, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )
def __Combine_Dictionaries__( dir_folder_mtx_10x_input, name_dict ) :
    """ # 2022-03-06 00:06:23 
    combined dictionaries processed from individual files
    """
    if os.path.exists( f'{dir_folder_mtx_10x_input}{name_dict}.tsv.gz' ) :
        ''' if an output file already exists, read the file and return the combined dictionary '''
        dict_combined = pd.read_csv( f'{dir_folder_mtx_10x_input}{name_dict}.tsv.gz', sep = '\t', header = None, index_col = 0 ).iloc[ :, 0 ].to_dict( )
    else :
        ''' combine summarized results '''
        l_dir_file = glob.glob( f"{dir_folder_mtx_10x_input}{name_dict}.*" )
        counter = collections.Counter( pd.read_csv( l_dir_file[ 0 ], sep = '\t', header = None, index_col = 0 ).iloc[ :, 0 ].to_dict( ) ) # initialize counter object with the dictionary from the first file
        for dir_file in l_dir_file[ 1 : ] :
            counter = counter + collections.Counter( pd.read_csv( dir_file, sep = '\t', header = None, index_col = 0 ).iloc[ :, 0 ].to_dict( ) ) # update counter object using the dictionary from each file
        dict_combined = dict( counter ) # retrieve a combined dictionary
        '''remove temporary files '''
        for dir_file in l_dir_file :
            os.remove( dir_file )
        ''' save dictionary as a file '''
        pd.Series( dict_combined ).to_csv( f'{dir_folder_mtx_10x_input}{name_dict}.tsv.gz', sep = '\t', header = None )
    return dict_combined # returns a combined dictionary
    
def __MTX_10X_Summarize_Counts__summarize_counts_for_each_mtx_10x__( dir_file_input, dir_folder_mtx_10x_input ) :
    '''
    internal function for MTX_10X_Filter
    # 2022-02-17 21:26:32   
    '''
    ''' survey the metrics '''
    ''' for each split mtx file, count number of umi and n_feature for each cells or the number of cells for each feature '''
    ''' initialize the dictionaries that will be handled by the current function '''
    dict_id_column_to_count = dict( )
    dict_id_column_to_n_features = dict( )
    dict_id_row_to_count = dict( )
    dict_id_row_to_n_cells = dict( )
    dict_id_row_to_log_transformed_count = dict( )
    
    for dir_file_input_mtx in pd.read_csv( dir_file_input, sep = '\t', header = None ).values.ravel( ) :
        with gzip.open( dir_file_input_mtx, 'rb' ) as file :
            ''' read the first line '''
            line = file.readline( ).decode( ) 
            ''' if the first line of the file contains a comment line, read all comment lines and a description line following the comments. '''
            if line[ 0 ] == '%' :
                # read comment and the description line
                while True :
                    if line[ 0 ] != '%' :
                        break
                    line = file.readline( ).decode( ) # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
                line = file.readline( ).decode( ) # read the next line
            ''' process entries'''
            while True :
                if len( line ) == 0 :
                    break
                ''' parse a record, and update metrics '''
                id_row, id_column, int_value = tuple( int( e ) for e in line.strip( ).split( ) ) # parse a record of a matrix-market format file
                ''' 1-based > 0-based coordinates '''
                id_row -= 1
                id_column -= 1
                ''' update umi count for each cell '''
                if id_column not in dict_id_column_to_count :
                    dict_id_column_to_count[ id_column ] = 0
                dict_id_column_to_count[ id_column ] += int_value
                ''' update n_features for each cell '''
                if id_column not in dict_id_column_to_n_features :
                    dict_id_column_to_n_features[ id_column ] = 0
                dict_id_column_to_n_features[ id_column ] += 1
                ''' update umi count for each feature '''
                if id_row not in dict_id_row_to_count :
                    dict_id_row_to_count[ id_row ] = 0
                dict_id_row_to_count[ id_row ] += int_value
                ''' update n_cells for each feature '''
                if id_row not in dict_id_row_to_n_cells :
                    dict_id_row_to_n_cells[ id_row ] = 0
                dict_id_row_to_n_cells[ id_row ] += 1
                ''' update log transformed counts, calculated by 'X_new = log_10(X_old + 1)', for each feature '''
                if id_row not in dict_id_row_to_log_transformed_count :
                    dict_id_row_to_log_transformed_count[ id_row ] = 0
                dict_id_row_to_log_transformed_count[ id_row ] += math.log10( int_value + 1 )
                line = file.readline( ).decode( ) # binary > uncompressed string # read the next line
    
    # save collected count as tsv files
    str_uuid_process = UUID( ) # retrieve uuid of the current process
    pd.Series( dict_id_column_to_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_column_to_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_column_to_n_features ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_column_to_n_features.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_n_cells ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_n_cells.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_log_transformed_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_log_transformed_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
def MTX_10X_Summarize_Counts( dir_folder_mtx_10x_input, verbose = False, int_num_threads = 15, flag_split_mtx = True, int_max_num_entries_for_chunk = 10000000 ) :
    """ # 2022-02-23 22:54:35 
    Summarize 
    (1) UMI and Feature counts for cells, 
    (2) UMI and Cell counts for features, and
    (3) log10-transformed values of UMI counts (X_new = log_10(X_old + 1)) for features
    and save these metrics as TSV files
    
    Returns:
    dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count
    """

    ''' the name of the dictionaries handled by this function '''
    l_name_dict = [ 'dict_id_column_to_count', 'dict_id_column_to_n_features', 'dict_id_row_to_count', 'dict_id_row_to_n_cells', 'dict_id_row_to_log_transformed_count' ]
    
    ''' handle inputs '''
    if dir_folder_mtx_10x_input[ -1 ] != '/' :
        dir_folder_mtx_10x_input += '/'

    # define flag and check whether the flag exists
    dir_file_flag = f"{dir_folder_mtx_10x_input}counts_summarized.flag"
    if not os.path.exists( dir_file_flag ) :
        # define input file directories
        dir_file_input_bc = f'{dir_folder_mtx_10x_input}barcodes.tsv.gz'
        dir_file_input_feature = f'{dir_folder_mtx_10x_input}features.tsv.gz'
        dir_file_input_mtx = f'{dir_folder_mtx_10x_input}matrix.mtx.gz'

        # check whether all required files are present
        if sum( list( not os.path.exists( dir_folder ) for dir_folder in [ dir_file_input_bc, dir_file_input_feature, dir_file_input_mtx ] ) ) :
            if verbose :
                print( f'required file(s) is not present at {dir_folder_mtx_10x}' )

        ''' split input mtx file into multiple files '''
        l_dir_file_mtx_10x = MTX_10X_Split( dir_folder_mtx_10x_input, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )

        ''' summarize each split mtx file '''
        Multiprocessing( l_dir_file_mtx_10x, __MTX_10X_Summarize_Counts__summarize_counts_for_each_mtx_10x__, n_threads = int_num_threads, global_arguments = [ dir_folder_mtx_10x_input ] )

        ''' combine summarized results '''
        dict_dict = dict( )
        for name_dict in l_name_dict :
            dict_dict[ name_dict ] = __Combine_Dictionaries__( dir_folder_mtx_10x_input, name_dict )
        # write the flag
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed at ' + TIME_GET_timestamp( True ) )
    else :
        ''' read summarized results '''
        dict_dict = dict( )
        for name_dict in l_name_dict :
            dict_dict[ name_dict ] = pd.read_csv( f'{dir_folder_mtx_10x_input}{name_dict}.tsv.gz', sep = '\t', header = None, index_col = 0 ).iloc[ :, 0 ].to_dict( )
                          
    # return summarized metrics
    return tuple( dict_dict[ name_dict ] for name_dict in l_name_dict )
def MTX_10X_Retrieve_number_of_rows_columns_and_entries( dir_folder_mtx_10x_input ) :
    """ # 2022-03-05 19:58:32 
    Retrieve the number of rows, columns, and entries from the matrix with the matrix market format 
    
    Returns:
    int_num_rows, int_num_columns, int_num_entries
    """
    ''' handle inputs '''
    if dir_folder_mtx_10x_input[ -1 ] != '/' :
        dir_folder_mtx_10x_input += '/'

    # define input file directories
    dir_file_input_mtx = f'{dir_folder_mtx_10x_input}matrix.mtx.gz'
    
    # check whether all required files are present
    if sum( list( not os.path.exists( dir_folder ) for dir_folder in [ dir_file_input_mtx ] ) ) :
        if verbose :
            print( f'required file(s) is not present at {dir_folder_mtx_10x}' )
    
    # read the input matrix
    with gzip.open( dir_file_input_mtx, 'rb' ) as file :
        ''' read the first line '''
        line = file.readline( ).decode( ) 
        ''' if the first line of the file contains a comment line, read all comment lines and a description line following the comments. '''
        if line[ 0 ] == '%' :
            # read comment and the description line
            while True :
                if line[ 0 ] != '%' :
                    break
                line = file.readline( ).decode( ) # read the next line
            # process the description line
            int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
        else :
            ''' the first line does not contain a comment, assumes it contains a description line '''
            int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
    return int_num_rows, int_num_columns, int_num_entries

dict_id_column_to_count, dict_id_row_to_avg_count, dict_id_row_to_avg_log_transformed_count, dict_id_row_to_avg_normalized_count, dict_id_row_to_avg_log_transformed_normalized_count = dict( ), dict( ), dict( ), dict( ), dict( ) # global variables # total UMI counts for each cell, average feature counts for each feature
def __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__first_pass__( dir_file_input, dir_folder_mtx_10x_input, int_target_sum ) :
    ''' # 2022-03-06 01:21:07 
    internal function for MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr
    '''
    global dict_id_column_to_count, dict_id_row_to_avg_count, dict_id_row_to_avg_log_transformed_count # use data in read-only global variables 
    ''' initialize the dictionaries that will be handled by the current function '''
    dict_id_row_to_deviation_from_mean_count = dict( )
    dict_id_row_to_deviation_from_mean_log_transformed_count = dict( )
    dict_id_row_to_normalized_count = dict( )
    dict_id_row_to_log_transformed_normalized_count = dict( )
    
    for dir_file_input_mtx in pd.read_csv( dir_file_input, sep = '\t', header = None ).values.ravel( ) :
        with gzip.open( dir_file_input_mtx, 'rb' ) as file :
            ''' read the first line '''
            line = file.readline( ).decode( ) 
            ''' if the first line of the file contains a comment line, read all comment lines and a description line following the comments. '''
            if line[ 0 ] == '%' :
                # read comment and the description line
                while True :
                    if line[ 0 ] != '%' :
                        break
                    line = file.readline( ).decode( ) # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
                line = file.readline( ).decode( ) # read the next line
            ''' process entries'''
            while True :
                if len( line ) == 0 :
                    break
                ''' parse a record, and update metrics '''
                id_row, id_column, int_value = tuple( int( e ) for e in line.strip( ).split( ) ) # parse a record of a matrix-market format file
                ''' 1-based > 0-based coordinates '''
                id_row -= 1
                id_column -= 1
                
                ''' update deviation from mean umi count for count of each feature '''
                if id_row not in dict_id_row_to_deviation_from_mean_count :
                    dict_id_row_to_deviation_from_mean_count[ id_row ] = 0
                dict_id_row_to_deviation_from_mean_count[ id_row ] += ( int_value - dict_id_row_to_avg_count[ id_row ] ) ** 2
                ''' update deviation from mean log transformed umi count for log_transformed count of each feature '''
                if id_row not in dict_id_row_to_deviation_from_mean_log_transformed_count :
                    dict_id_row_to_deviation_from_mean_log_transformed_count[ id_row ] = 0
                dict_id_row_to_deviation_from_mean_log_transformed_count[ id_row ] += ( math.log10( int_value + 1 ) - dict_id_row_to_avg_log_transformed_count[ id_row ] ) ** 2
                ''' calculate normalized target sum '''
                int_value_normalized = int_value / dict_id_column_to_count[ id_column ] * int_target_sum 
                ''' update normalized counts, calculated by 'X_new = X_old / total_umi * int_target_sum', for each feature '''
                if id_row not in dict_id_row_to_normalized_count :
                    dict_id_row_to_normalized_count[ id_row ] = 0
                dict_id_row_to_normalized_count[ id_row ] += int_value_normalized
                ''' update log transformed normalized counts, calculated by 'X_new = log_10(X_old / total_umi * int_target_sum + 1)', for each feature '''
                if id_row not in dict_id_row_to_log_transformed_normalized_count :
                    dict_id_row_to_log_transformed_normalized_count[ id_row ] = 0
                dict_id_row_to_log_transformed_normalized_count[ id_row ] += math.log10( int_value_normalized + 1 ) 
                
                line = file.readline( ).decode( ) # binary > uncompressed string # read the next line
    
    # save collected count as tsv files
    str_uuid_process = UUID( ) # retrieve uuid of the current process
    pd.Series( dict_id_row_to_deviation_from_mean_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_deviation_from_mean_log_transformed_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_log_transformed_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_normalized_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_normalized_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_log_transformed_normalized_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_log_transformed_normalized_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
def __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__second_pass__( dir_file_input, dir_folder_mtx_10x_input, int_target_sum ) :
    ''' # 2022-03-06 01:21:14 
    internal function for MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr
    '''
    global dict_id_column_to_count, dict_id_row_to_avg_normalized_count, dict_id_row_to_avg_log_transformed_normalized_count # use data in read-only global variables 
    ''' initialize the dictionaries that will be handled by the current function '''
    dict_id_row_to_deviation_from_mean_normalized_count = dict( )
    dict_id_row_to_deviation_from_mean_log_transformed_normalized_count = dict( )
    
    for dir_file_input_mtx in pd.read_csv( dir_file_input, sep = '\t', header = None ).values.ravel( ) :
        with gzip.open( dir_file_input_mtx, 'rb' ) as file :
            ''' read the first line '''
            line = file.readline( ).decode( ) 
            ''' if the first line of the file contains a comment line, read all comment lines and a description line following the comments. '''
            if line[ 0 ] == '%' :
                # read comment and the description line
                while True :
                    if line[ 0 ] != '%' :
                        break
                    line = file.readline( ).decode( ) # read the next line
                # process the description line
                int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
                line = file.readline( ).decode( ) # read the next line
            ''' process entries'''
            while True :
                if len( line ) == 0 :
                    break
                ''' parse a record, and update metrics '''
                id_row, id_column, int_value = tuple( int( e ) for e in line.strip( ).split( ) ) # parse a record of a matrix-market format file
                ''' 1-based > 0-based coordinates '''
                id_row -= 1
                id_column -= 1
                
                ''' calculate normalized target sum '''
                int_value_normalized = int_value / dict_id_column_to_count[ id_column ] * int_target_sum 
                ''' update deviation from mean normalized umi counts, calculated by 'X_new = X_old / total_umi * int_target_sum', for each feature '''
                if id_row not in dict_id_row_to_deviation_from_mean_normalized_count :
                    dict_id_row_to_deviation_from_mean_normalized_count[ id_row ] = 0
                dict_id_row_to_deviation_from_mean_normalized_count[ id_row ] += ( int_value_normalized - dict_id_row_to_avg_normalized_count[ id_row ] ) ** 2
                ''' update deviation from mean log transformed normalized umi counts, calculated by 'X_new = log_10(X_old / total_umi * int_target_sum + 1)', for each feature '''
                if id_row not in dict_id_row_to_deviation_from_mean_log_transformed_normalized_count :
                    dict_id_row_to_deviation_from_mean_log_transformed_normalized_count[ id_row ] = 0
                dict_id_row_to_deviation_from_mean_log_transformed_normalized_count[ id_row ] += ( math.log10( int_value_normalized + 1 ) - dict_id_row_to_avg_log_transformed_normalized_count[ id_row ] ) ** 2
                
                line = file.readline( ).decode( ) # binary > uncompressed string # read the next line
    
    # save collected count as tsv files
    str_uuid_process = UUID( ) # retrieve uuid of the current process
    pd.Series( dict_id_row_to_deviation_from_mean_normalized_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_normalized_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
    pd.Series( dict_id_row_to_deviation_from_mean_log_transformed_normalized_count ).to_csv( f'{dir_folder_mtx_10x_input}dict_id_row_to_deviation_from_mean_log_transformed_normalized_count.{str_uuid_process}.tsv.gz', sep = '\t', header = None )
def MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr( dir_folder_mtx_10x_input, int_target_sum = 10000, verbose = False, int_num_threads = 15, flag_split_mtx = True, int_max_num_entries_for_chunk = 10000000 ) :
    """ # 2022-02-23 22:54:35 
    Calculate average log transformed normalized expr
    (1) UMI and Feature counts for cells, and
    (2) Cell counts for features,
    and save these metrics as TSV files
    
    Arguments:
    'int_target_sum' : the target count for the total UMI count for each cell. The counts will normalized to meet the target sum.
    
    Returns:
    dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count
    """

    ''' handle inputs '''
    if dir_folder_mtx_10x_input[ -1 ] != '/' :
        dir_folder_mtx_10x_input += '/'

    # define flag and check whether the flag exists
    dir_file_flag = f"{dir_folder_mtx_10x_input}avg_expr_normalized_summarized.int_target_sum__{int_target_sum}.flag"
    if not os.path.exists( dir_file_flag ) :
        # define input file directories
        dir_file_input_bc = f'{dir_folder_mtx_10x_input}barcodes.tsv.gz'
        dir_file_input_feature = f'{dir_folder_mtx_10x_input}features.tsv.gz'
        dir_file_input_mtx = f'{dir_folder_mtx_10x_input}matrix.mtx.gz'

        # check whether all required files are present
        if sum( list( not os.path.exists( dir_folder ) for dir_folder in [ dir_file_input_bc, dir_file_input_feature, dir_file_input_mtx ] ) ) :
            if verbose :
                print( f'required file(s) is not present at {dir_folder_mtx_10x}' )

        ''' split input mtx file into multiple files '''
        l_dir_file_mtx_10x = MTX_10X_Split( dir_folder_mtx_10x_input, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )

        ''' retrieve number of cells, features, and entries from the matrix file '''
        int_num_cells, int_num_features, int_num_entries = MTX_10X_Retrieve_number_of_rows_columns_and_entries( dir_folder_mtx_10x_input )
        
        ''' summarizes counts '''
        global dict_id_column_to_count, dict_id_row_to_avg_count, dict_id_row_to_avg_log_transformed_count, dict_id_row_to_avg_normalized_count, dict_id_row_to_avg_log_transformed_normalized_count # use global variable
        dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count = MTX_10X_Summarize_Counts( dir_folder_mtx_10x_input, verbose = verbose, int_num_threads = int_num_threads, flag_split_mtx = flag_split_mtx, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )

        """ first pass """
        # calculate mean counts
        dict_id_row_to_avg_count = ( pd.Series( dict_id_row_to_count ) / int_num_cells ).to_dict( ) # calculate average expression of each feature
        dict_id_row_to_avg_log_transformed_count = ( pd.Series( dict_id_row_to_log_transformed_count ) / int_num_cells ).to_dict( ) # calculate average log-transformed expression of each feature
        
        ''' calculated average log2 transformed normalized expr for each split mtx file '''
        Multiprocessing( l_dir_file_mtx_10x, __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__first_pass__, n_threads = int_num_threads, global_arguments = [ dir_folder_mtx_10x_input, int_target_sum ] )

        l_name_dict_first_pass = [ 'dict_id_row_to_deviation_from_mean_count', 'dict_id_row_to_deviation_from_mean_log_transformed_count', 'dict_id_row_to_normalized_count', 'dict_id_row_to_log_transformed_normalized_count' ]
        
        ''' combine summarized results '''
        dict_dict = dict( )
        for name_dict in l_name_dict_first_pass :
            dict_dict[ name_dict ] = __Combine_Dictionaries__( dir_folder_mtx_10x_input, name_dict )
            
        """ second pass """
        # calculate mean counts
        dict_id_row_to_avg_normalized_count = ( pd.Series( dict_dict[ 'dict_id_row_to_normalized_count' ] ) / int_num_cells ).to_dict( ) # calculate average expression of each feature
        dict_id_row_to_avg_log_transformed_normalized_count = ( pd.Series( dict_dict[ 'dict_id_row_to_log_transformed_normalized_count' ] ) / int_num_cells ).to_dict( ) # calculate average log-transformed expression of each feature
        
        ''' calculated average log2 transformed normalized expr for each split mtx file '''
        Multiprocessing( l_dir_file_mtx_10x, __MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr__second_pass__, n_threads = int_num_threads, global_arguments = [ dir_folder_mtx_10x_input, int_target_sum ] )

        l_name_dict_second_pass = [ 'dict_id_row_to_deviation_from_mean_normalized_count', 'dict_id_row_to_deviation_from_mean_log_transformed_normalized_count' ]
        
        ''' combine summarized results '''
        for name_dict in l_name_dict_second_pass :
            dict_dict[ name_dict ] = __Combine_Dictionaries__( dir_folder_mtx_10x_input, name_dict )
            
        ''' compose a dataframe containing the summary about the features '''
        df_summary = pd.DataFrame( { 
            'n_cells' : pd.Series( dict_id_row_to_n_cells ),
            'variance_of_count' : pd.Series( dict_dict[ 'dict_id_row_to_deviation_from_mean_count' ] ) / ( int_num_cells - 1 ),
            'variance_of_log_transformed_count' : pd.Series( dict_dict[ 'dict_id_row_to_deviation_from_mean_log_transformed_count' ] ) / ( int_num_cells - 1 ),
            'variance_of_normalized_count' : pd.Series( dict_dict[ 'dict_id_row_to_deviation_from_mean_normalized_count' ] ) / ( int_num_cells - 1 ),
            'variance_of_log_transformed_normalized_count' : pd.Series( dict_dict[ 'dict_id_row_to_deviation_from_mean_log_transformed_normalized_count' ] ) / ( int_num_cells - 1 ),
            'mean_count' : pd.Series( dict_id_row_to_avg_count ),
            'mean_log_transformed_count' : pd.Series( dict_id_row_to_avg_log_transformed_count ),
            'mean_normalized_count' : pd.Series( dict_id_row_to_avg_normalized_count ),
            'mean_log_transformed_normalized_count' : pd.Series( dict_id_row_to_avg_log_transformed_normalized_count ),
        } )
        # read a dataframe containing features
        df_feature = pd.read_csv( dir_file_input_feature, sep = '\t', header = None )
        df_feature.columns = [ 'id_feature', 'feature', '10X_type' ]
        
        df_summary = df_summary.join( df_feature, how = 'left' ) # add df_feature to the df_summary
        df_summary.index.name = 'id_row' 
        df_summary.reset_index( drop = False, inplace = True ) # retrieve id_row as a column
        df_summary.to_csv( f'{dir_folder_mtx_10x_input}statistical_summary_of_features.int_target_sum__{int_target_sum}.tsv.gz', sep = '\t', index = False ) # save statistical summary as a text file
        
        # write the flag
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed at ' + TIME_GET_timestamp( True ) )
    else :
        ''' if 'MTX_10X_Calculate_Average_Log10_Transformed_Normalized_Expr' function has been already run on the current folder, read the previously saved result, and return the summary dataframe '''
        df_summary = pd.read_csv( f'{dir_folder_mtx_10x_input}statistical_summary_of_features.int_target_sum__{int_target_sum}.tsv.gz', sep = '\t' ) # save statistical summary as a text file
    return df_summary
dict_id_column_previous_to_id_column_current, dict_id_row_previous_to_id_row_current = dict( ), dict( )
def __MTX_10X_Filter__filter_mtx_10x__( dir_file_input, dir_folder_mtx_10x_output ) :
    """ # 2022-02-22 02:06:03 
    __MTX_10X_Filter__filter_mtx_10x__
    
    Returns:
    int_n_entries = total number of entries written by the current process after filtering
    """
    int_n_entries = 0 # total number of entries written by the current process after filtering
#     dict_id_column_previous_to_id_column_current = PICKLE_Read( f'{dir_folder_mtx_10x_output}dict_id_column_previous_to_id_column_current.pickle' ) # retrieve id_feature to index_feature mapping 
#     dict_id_row_previous_to_id_row_current = PICKLE_Read( f'{dir_folder_mtx_10x_output}dict_id_row_previous_to_id_row_current.pickle' ) # retrieve id_feature to index_feature mapping 
    """ write a filtered matrix.mtx.gz for each split mtx file """
    for dir_file_mtx_10x, index_mtx_10x in pd.read_csv( dir_file_input, sep = '\t' ).values :
        # directly write matrix.mtx.gz file without using an external dependency
        with gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.{index_mtx_10x}.gz", 'wb' ) as newfile :
            with gzip.open( dir_file_mtx_10x, 'rb' ) as file : 
                ''' read the first line '''
                line = file.readline( ).decode( ) 
                ''' if the first line of the file contains a comment line, read all comment lines and a description line following the comments. '''
                if line[ 0 ] == '%' :
                    # read comment and the description line
                    while True :
                        if line[ 0 ] != '%' :
                            break
                        line = file.readline( ).decode( ) # read the next line
                    # process the description line
                    int_num_rows, int_num_columns, int_num_entries = tuple( int( e ) for e in line.strip( ).split( ) ) # retrieve the number of rows, number of columns and number of entries
                    line = file.readline( ).decode( ) # read the next line
                ''' process entries'''
                while True :
                    if len( line ) == 0 :
                        break
                    id_row, id_column, int_value = tuple( map( int, line.strip( ).split( ) ) ) # parse each entry of the current matrix 
                    ''' 1-based > 0-based coordinates '''
                    id_row -= 1
                    id_column -= 1
                    ''' write a record to the new matrix file only when both id_row and id_column belongs to filtered id_rows and id_columns '''
                    if id_row in dict_id_row_previous_to_id_row_current and id_column in dict_id_column_previous_to_id_column_current :
                        newfile.write( ( ' '.join( tuple( map( str, [ dict_id_row_previous_to_id_row_current[ id_row ] + 1, dict_id_column_previous_to_id_column_current[ id_column ] + 1, int_value ] ) ) ) + '\n' ).encode( ) ) # map id_row and id_column of the previous matrix to those of the filtered matrix (new matrix) # 0-based > 1-based coordinates
                        int_n_entries += 1 # update the total number of entries written by the current process
                    line = file.readline( ).decode( ) # read the next line
    return int_n_entries # returns the total number of entries written by the current process
def MTX_10X_Filter( dir_folder_mtx_10x_input, dir_folder_mtx_10x_output, min_counts = None, min_features = None, min_cells = None, l_features = None, l_cells = None, verbose = False, function_for_adjusting_thresholds = None, int_num_threads = 15, flag_split_mtx = True, int_max_num_entries_for_chunk = 10000000 ) :
    ''' # 2022-02-22 01:39:45  hyunsu-an
    read 10x count matrix and filter matrix based on several thresholds
    'dir_folder_mtx_10x_input' : a folder containing files for the input 10x count matrix
    'dir_folder_mtx_10x_output' : a folder containing files for the input 10x count matrix

    Only the threshold arguments for either cells ( 'min_counts', 'min_features' ) or features ( 'min_cells' ) can be given at a time.

    'min_counts' : the minimum number of total counts for a cell to be included in the output matrix
    'min_features' : the minimum number of features for a cell to be included in the output matrix
    'min_cells' : the minimum number of cells for a feature to be included in the output matrix
    'l_features' : a list of features (values in the first column of 'features.tsv.gz') to include. All other features will be excluded from the output matrix. (default: None) If None is given, include all features in the output matrix.
    'l_cells' : a list of cells (values in the first column of 'barcodes.tsv.gz') to include. All other cells will be excluded from the output matrix. (default: None) If None is given, include all cells in the output matrix.
    'int_num_threads' : when 'int_num_threads' is 1, does not use the multiprocessing  module for parallel processing
    'function_for_adjusting_thresholds' : a function for adjusting thresholds based on the summarized metrics. Useful when the exact threshold for removing empty droplets are variable across the samples. the function should receive arguments and return values in the following structure: 
                                        min_counts_new, min_features_new, min_cells_new = function_for_adjusting_thresholds( dir_folder_mtx_10x_output, min_counts, min_features, min_cells )
    '''

    ''' handle inputs '''
    if dir_folder_mtx_10x_input[ -1 ] != '/' :
        dir_folder_mtx_10x_input += '/'
    if dir_folder_mtx_10x_output[ -1 ] != '/' :
        dir_folder_mtx_10x_output += '/'
    if ( ( min_counts is not None ) or ( min_features is not None ) ) and ( min_cells is not None ) : # check whether thresholds for both cells and features were given (thresdholds for either cells or features can be given at a time)
        if verbose :
            print( '[MTX_10X_Filter] (error) no threshold is given or more thresholds for both cells and features are given. (Thresdholds for either cells or features can be given at a time.)' )
        return -1
    # create an output folder
    os.makedirs( dir_folder_mtx_10x_output, exist_ok = True )

    # define input file directories
    dir_file_input_bc = f'{dir_folder_mtx_10x_input}barcodes.tsv.gz'
    dir_file_input_feature = f'{dir_folder_mtx_10x_input}features.tsv.gz'
    dir_file_input_mtx = f'{dir_folder_mtx_10x_input}matrix.mtx.gz'

    # check whether all required files are present
    if sum( list( not os.path.exists( dir_folder ) for dir_folder in [ dir_file_input_bc, dir_file_input_feature, dir_file_input_mtx ] ) ) :
        if verbose :
            print( f'required file(s) is not present at {dir_folder_mtx_10x}' )

    ''' read barcode and feature information '''
    df_bc = pd.read_csv( dir_file_input_bc, sep = '\t', header = None )
    df_bc.columns = [ 'barcode' ]
    df_feature = pd.read_csv( dir_file_input_feature, sep = '\t', header = None )
    df_feature.columns = [ 'id_feature', 'feature', '10X_type' ]

    ''' split input mtx file into multiple files '''
    l_dir_file_mtx_10x = MTX_10X_Split( dir_folder_mtx_10x_input, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )
    
    ''' summarizes counts '''
    dict_id_column_to_count, dict_id_column_to_n_features, dict_id_row_to_count, dict_id_row_to_n_cells, dict_id_row_to_log_transformed_count = MTX_10X_Summarize_Counts( dir_folder_mtx_10x_input, verbose = verbose, int_num_threads = int_num_threads, flag_split_mtx = flag_split_mtx, int_max_num_entries_for_chunk = int_max_num_entries_for_chunk )
    
    ''' adjust thresholds based on the summarized metrices (if a function has been given) '''
    if function_for_adjusting_thresholds is not None :
        min_counts, min_features, min_cells = function_for_adjusting_thresholds( dir_folder_mtx_10x_output, min_counts, min_features, min_cells )
    
    ''' filter row or column that do not satisfy the given thresholds '''
    if min_counts is not None :
        dict_id_column_to_count = dict( ( k, dict_id_column_to_count[ k ] ) for k in dict_id_column_to_count if dict_id_column_to_count[ k ] >= min_counts ) 
    if min_features is not None :
        dict_id_column_to_n_features = dict( ( k, dict_id_column_to_n_features[ k ] ) for k in dict_id_column_to_n_features if dict_id_column_to_n_features[ k ] >= min_features )
    if min_cells is not None :
        dict_id_row_to_n_cells = dict( ( k, dict_id_row_to_n_cells[ k ] ) for k in dict_id_row_to_n_cells if dict_id_row_to_n_cells[ k ] >= min_cells )

    ''' retrieve id_row and id_column that satisfy the given thresholds '''    
    set_id_column = set( dict_id_column_to_count ).intersection( set( dict_id_column_to_n_features ) )
    set_id_row = set( dict_id_row_to_n_cells )
    
    ''' exclude cells and features not present in the input lists (if the lists were given)  '''
    if l_cells is not None :        
        dict_barcode_to_id_column = dict( ( barcode, id_column ) for id_column, barcode in enumerate( df_bc.barcode.values ) )
        set_id_column = set_id_column.intersection( set( dict_barcode_to_id_column[ barcode ] for barcode in set( l_cells ) if barcode in dict_barcode_to_id_column ) )
        del dict_barcode_to_id_column
    if l_features is not None :
        dict_id_feature_to_id_row = dict( ( id_feature, id_row ) for id_row, id_feature in enumerate( df_feature.id_feature.values ) )
        set_id_row = set_id_row.intersection( set( dict_id_feature_to_id_row[ id_feature ] for id_feature in set( l_features ) if id_feature in dict_id_feature_to_id_row ) )
        del dict_id_feature_to_id_row

    ''' report the number of cells or features that will be filtered out '''
    if verbose :
        int_n_bc_filtered = len( df_bc ) - len( set_id_column )
        if int_n_bc_filtered > 0 :
            print( f"{int_n_bc_filtered}/{len( df_bc )} barcodes will be filtered out" )
        int_n_feature_filtered = len( df_feature ) - len( set_id_row )
        if int_n_feature_filtered > 0 :
            print( f"{int_n_feature_filtered}/{len( df_feature )} features will be filtered out" )

    """ retrieve a mapping between previous id_column to current id_column """
    global dict_id_column_previous_to_id_column_current, dict_id_row_previous_to_id_row_current # use global variables for multiprocessing
    df_bc = df_bc.loc[ list( set_id_column ) ]
    df_bc.index.name = 'id_column_previous'
    df_bc.reset_index( drop = False, inplace = True )
    df_bc[ 'id_column_current' ] = np.arange( len( df_bc ) )
    dict_id_column_previous_to_id_column_current = df_bc.set_index( 'id_column_previous' ).id_column_current.to_dict( ) 
    PICKLE_Write( f'{dir_folder_mtx_10x_output}dict_id_column_previous_to_id_column_current.pickle', dict_id_column_previous_to_id_column_current ) # save id_feature to index_feature mapping 
    """ retrieve a mapping between previous id_row to current id_row """
    df_feature = df_feature.loc[ list( set_id_row ) ]
    df_feature.index.name = 'id_row_previous'
    df_feature.reset_index( drop = False, inplace = True )
    df_feature[ 'id_row_current' ] = np.arange( len( df_feature ) )
    dict_id_row_previous_to_id_row_current = df_feature.set_index( 'id_row_previous' ).id_row_current.to_dict( ) 
    PICKLE_Write( f'{dir_folder_mtx_10x_output}dict_id_row_previous_to_id_row_current.pickle', dict_id_row_previous_to_id_row_current ) # save id_feature to index_feature mapping 

    ''' save barcode file '''
    df_bc.to_csv( f"{dir_folder_mtx_10x_output}barcodes.tsv.gz", columns = [ 'barcode' ], sep = '\t', index = False, header = False ) 

    ''' save feature file '''
    df_feature[ [ 'id_feature', 'feature', '10X_type' ] ].to_csv( f"{dir_folder_mtx_10x_output}features.tsv.gz", sep = '\t', index = False, header = False ) # save as a file

    """ write a filtered matrix.mtx.gz for each split mtx file using multiple processes and retrieve the total number of entries written by each process """
    # compose inputs for multiprocessing
    df_input = pd.DataFrame( { 'dir_file_mtx_10x' : l_dir_file_mtx_10x, 'index_mtx_10x' : np.arange( len( l_dir_file_mtx_10x ) ) } )
    l_int_n_entries = Multiprocessing( df_input, __MTX_10X_Filter__filter_mtx_10x__, int_num_threads, global_arguments = [ dir_folder_mtx_10x_output ] ) 
    # retrieve the total number of entries
    int_total_n_entries = sum( l_int_n_entries )
    
    """ combine parts and add the MTX file header to compose a combined mtx file """
    df_file = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.*.gz" )
    df_file.wildcard_0 = df_file.wildcard_0.astype( int )
    df_file.sort_values( 'wildcard_0', inplace = True )
    OS_FILE_Combine_Files_in_order( df_file.dir.values, f"{dir_folder_mtx_10x_output}matrix.mtx.gz", delete_input_files = not flag_split_mtx, header = f"%%MatrixMarket matrix coordinate integer general\n%\n{len( dict_id_row_previous_to_id_row_current )} {len( dict_id_column_previous_to_id_column_current )} {int_total_n_entries}\n" ) # combine the output mtx files in the order # does not delete temporary files if 'flag_split_mtx' is True
    
    # write a flag indicating that the current output directory contains split mtx files
    with open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz.split.flag", 'w' ) as file :
        file.write( 'completed' )