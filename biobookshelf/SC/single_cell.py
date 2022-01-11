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
    
def MTX_10X_Filter( dir_folder_mtx_10x_input, dir_folder_mtx_10x_output, min_counts = None, min_features = None, min_cells = None, verbose = False ) :
    ''' # 2022-01-08 20:25:37 hyunsu-an
    read 10x count matrix and filter matrix based on several thresholds
    'dir_folder_mtx_10x_input' : a folder containing files for the input 10x count matrix
    'dir_folder_mtx_10x_output' : a folder containing files for the input 10x count matrix

    Only the threshold arguments for either cells ( 'min_counts', 'min_features' ) or features ( 'min_cells' ) can be given at a time.

    'min_counts' : the minimum number of total counts for a cell to be included in the output matrix
    'min_features' : the minimum number of features for a cell to be included in the output matrix
    'min_cells' : the minimum number of cells for a feature to be included in the output matrix
    '''

    ''' handle inputs '''
    if dir_folder_mtx_10x_input[ -1 ] != '/' :
        dir_folder_mtx_10x_input += '/'
    if dir_folder_mtx_10x_output[ -1 ] != '/' :
        dir_folder_mtx_10x_output += '/'
    if not ( ( ( min_counts is not None ) or ( min_features is not None ) ) ^ ( min_cells is not None ) ) : # check whether a any threshold is given, or thresholds for both cells and features were given (thresdholds for either cells or features can be given at a time)
        if verbose :
            print( '[MTX_10X_Filter] (error) no threshold is given or more thresholds for both cells and features are given. (Thresdholds for either cells or features can be given at a time.)' )
        return -1
    #     pass

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

    ''' survey the metrics '''
    dict_id_column_to_count = dict( )
    dict_id_column_to_n_features = dict( )
    dict_id_row_to_n_cells = dict( )

    with gzip.open( dir_file_input_mtx, 'rb' ) as file :
        while True :
            line = file.readline( ).decode( ) # binary > uncompressed string
            if len( line ) == 0 :
                break
            ''' skip comments '''
            if line[ 0 ] == '%' :
                continue
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
            ''' update n_cells for each feature '''
            if id_row not in dict_id_row_to_n_cells :
                dict_id_row_to_n_cells[ id_row ] = 0
            dict_id_row_to_n_cells[ id_row ] += 1


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

    ''' report the number of cells or features that will be filtered out '''
    if verbose :
        int_n_bc_filtered = len( df_bc ) - len( set_id_column )
        if int_n_bc_filtered > 0 :
            print( f"{int_n_bc_filtered}/{len( df_bc )} barcodes will be filtered out" )
        int_n_feature_filtered = len( df_feature ) - len( set_id_row )
        if int_n_feature_filtered > 0 :
            print( f"{int_n_feature_filtered}/{len( df_feature )} features will be filtered out" )

    """ retrieve a mapping between previous id_column to current id_column """
    df_bc = df_bc.loc[ set_id_column ]
    df_bc.index.name = 'id_column_previous'
    df_bc.reset_index( drop = False, inplace = True )
    df_bc[ 'id_column_current' ] = np.arange( len( df_bc ) )
    dict_id_column_previous_to_id_column_current = df_bc.set_index( 'id_column_previous' ).id_column_current.to_dict( ) 
    """ retrieve a mapping between previous id_row to current id_row """
    df_feature = df_feature.loc[ set_id_row ]
    df_feature.index.name = 'id_row_previous'
    df_feature.reset_index( drop = False, inplace = True )
    df_feature[ 'id_row_current' ] = np.arange( len( df_feature ) )
    dict_id_row_previous_to_id_row_current = df_feature.set_index( 'id_row_previous' ).id_row_current.to_dict( ) 


    """ write a combined matrix.mtx.gz """
    # create an output folder
    os.makedirs( dir_folder_mtx_10x_output, exist_ok = True ) 
    # retrieve the total number of entries
    int_total_n_entries = pd.Series( dict( ( k, dict_id_column_to_n_features[ k ] ) for k in dict_id_column_to_n_features if k in set_id_column ) ).sum( ) # the total number of records after filtering can be retrieved by calculating the sum of values of the 'dict_id_column_to_n_features'
    # directly write matrix.mtx.gz file without using an external dependency
    with gzip.open( f"{dir_folder_mtx_10x_output}matrix.mtx.gz", 'wb' ) as newfile :
        newfile.write( ( f"%%MatrixMarket matrix coordinate integer general\n%\n{len( dict_id_row_previous_to_id_row_current )} {len( dict_id_column_previous_to_id_column_current )} {int_total_n_entries}\n" ).encode( ) ) # write matrix file header
        with gzip.open( dir_file_input_mtx, 'rb' ) as file : # retrieve a list of features
            while True :
                line = file.readline( ).decode( )
                if len( line ) == 0 :
                    break
                ''' skip comments '''
                if line[ 0 ] == '%' :
                    continue
                id_row, id_column, int_value = tuple( map( int, line.strip( ).split( ) ) ) # parse each entry of the current matrix 
                ''' 1-based > 0-based coordinates '''
                id_row -= 1
                id_column -= 1
                ''' write a record to the new matrix file only when both id_row and id_column belongs to filtered id_rows and id_columns '''
                if id_row in dict_id_row_previous_to_id_row_current and id_column in dict_id_column_previous_to_id_column_current :
                    newfile.write( ( ' '.join( tuple( map( str, [ dict_id_row_previous_to_id_row_current[ id_row ] + 1, dict_id_column_previous_to_id_column_current[ id_column ] + 1, int_value ] ) ) ) + '\n' ).encode( ) ) # map id_row and id_column of the previous matrix to those of the filtered matrix (new matrix) # 0-based > 1-based coordinates

    ''' save barcode file '''
    df_bc.to_csv( f"{dir_folder_mtx_10x_output}barcodes.tsv.gz", columns = [ 'barcode' ], sep = '\t', index = False, header = False ) 

    ''' save feature file '''
    df_feature[ [ 'id_feature', 'feature', '10X_type' ] ].to_csv( f"{dir_folder_mtx_10x_output}features.tsv.gz", sep = '\t', index = False, header = False ) # save as a file
    
def MTX_10X_Combine( dir_folder_output_mtx_10x, * l_dir_folder_input_mtx_10x ) :
    '''
    # 2021-12-18 13:13:03 
    Combine 10X count matrix files from the given list of folders and write combined output files to the given output folder 'dir_folder_output_mtx_10x'
    If there are no shared cells between matrix files, a low-memory mode will be used. The output files will be simply combined since no count summing operation is needed. Only feature matrix will be loaded and updated in the memory.
    'id_feature' should be unique across all features
    '''
    
    # create an output folder
    os.makedirs( dir_folder_output_mtx_10x, exist_ok = True ) 

    """ retrieve cell barcodes of all 10X matrices and check whether cell barcodes are not shared between matrices """
    dict_dir_folder_mtx_10x_to_int_n_barcodes = dict( ) # collect the number of barcodes for each matrix
    set_barcode = set( ) # update a set of unique barcodes
    for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x :
        arr_barcode = pd.read_csv( f'{dir_folder_mtx_10x}barcodes.tsv.gz', sep = '\t', header = None, squeeze = True ).values # retrieve a list of features
        set_barcode.update( arr_barcode ) # update a set of barcodes
        dict_dir_folder_mtx_10x_to_int_n_barcodes[ dir_folder_mtx_10x ] = len( arr_barcode ) # collect the number of barcodes 

    ''' check whether there are shared cell barcodes between matrices and set a flag for entering a low-memory mode '''
    flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs = len( set_barcode ) == sum( list( dict_dir_folder_mtx_10x_to_int_n_barcodes[ k ] for k in dict_dir_folder_mtx_10x_to_int_n_barcodes ) ) # update flag

    if flag_low_memory_mode_because_there_is_no_shared_cell_between_mtxs :
        """ low-memory mode """
        """ write a combined barcodes.tsv.gz """
        OS_FILE_Combine_Files_in_order( list( f"{dir_folder_mtx_10x}barcodes.tsv.gz" for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x ), f"{dir_folder_output_mtx_10x}barcodes.tsv.gz", overwrite_existing_file = True )

        ''' collect a set of unique features and a list of id_feature for each 10X matrix '''
        # assumes id_feature is unique for all features across 10X matrices
        dict_dir_folder_mtx_10x_to_arr_id_feature = dict( ) # collect a list of id_feature for each matrix folder
        set_t_feature = set( ) # update a set of feature tuples
        for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x :
            arr_feature = pd.read_csv( f'{dir_folder_mtx_10x}features.tsv.gz', sep = '\t', header = None ).values # retrieve a list of features
            set_t_feature.update( list( map( tuple, arr_feature ) ) ) # update a set of feature tuples
            dict_dir_folder_mtx_10x_to_arr_id_feature[ dir_folder_mtx_10x ] = arr_feature[ :, 0 ] # retrieve a list of id_feature

        """ write a combined features.tsv.gz """
        l_t_feature = list( set_t_feature ) # convert set to list
        with gzip.open( f"{dir_folder_output_mtx_10x}features.tsv.gz", 'wb' ) as newfile :
            for t_feature in l_t_feature :
                newfile.write( ( '\t'.join( t_feature ) + '\n' ).encode( ) )

        """ build a mapping of id_feature to index_feature, which will be consistent across datasets """
        dict_id_feature_to_index_feature = dict( ( t_feature[ 0 ], index_feature + 1 ) for index_feature, t_feature in enumerate( l_t_feature ) )

        ''' collect the number of entries for each 10X matrix '''
        dict_dir_folder_mtx_10x_to_int_n_entries = dict( ) # collect the number of entries for each matrix
        for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x :
            with gzip.open( f'{dir_folder_mtx_10x}matrix.mtx.gz', 'rb' ) as file : # retrieve a list of features
                file.readline( ), file.readline( )
                dict_dir_folder_mtx_10x_to_int_n_entries[ dir_folder_mtx_10x ] = int( file.readline( ).decode( ).strip( ).split( )[ 2 ] )
        int_total_n_entries = sum( list( dict_dir_folder_mtx_10x_to_int_n_entries[ k ] for k in dict_dir_folder_mtx_10x_to_int_n_entries ) ) # retrieve a total number of entries

        """ write a combined matrix.mtx.gz """
        # directly write matrix.mtx.gz file without using external dependency
        int_total_n_barcodes_of_previously_written_matrices = 0
        with gzip.open( f"{dir_folder_output_mtx_10x}matrix.mtx.gz", 'wb' ) as newfile :
            newfile.write( ( f"%%MatrixMarket matrix coordinate integer general\n%\n{len( l_t_feature )} {len( set_barcode )} {int_total_n_entries}\n" ).encode( ) ) # write matrix file header
            for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x :
                with gzip.open( f'{dir_folder_mtx_10x}matrix.mtx.gz', 'rb' ) as file : # retrieve a list of features
                    file.readline( ), file.readline( ), file.readline( ) # read three header lines
                    while True :
                        line = file.readline( ).decode( )
                        if len( line ) == 0 :
                            break
                        index_row, index_col, int_value = tuple( map( int, line.strip( ).split( ) ) ) # parse each entry of the current matrix 
                        newfile.write( ( ' '.join( tuple( map( str, [ dict_id_feature_to_index_feature[ dict_dir_folder_mtx_10x_to_arr_id_feature[ dir_folder_mtx_10x ][ index_row - 1 ] ], index_col + int_total_n_barcodes_of_previously_written_matrices, int_value ] ) ) ) + '\n' ).encode( ) ) # translate indices of the current matrix to that of the combined matrix
                int_total_n_barcodes_of_previously_written_matrices += dict_dir_folder_mtx_10x_to_int_n_barcodes[ dir_folder_mtx_10x ] # update the total number of barcodes of matrices that were completely written to the combined output matrix file
    else :
        ''' normal operation mode perfoming count merging operations '''
        l_df_mtx, l_df_feature = [ ], [ ]
        for dir_folder_mtx_10x in l_dir_folder_input_mtx_10x :
            df_mtx, df_feature = MTX_10X_Read( dir_folder_mtx_10x )
            l_df_mtx.append( df_mtx ), l_df_feature.append( df_feature )

        # combine mtx
        df_mtx = pd.concat( l_df_mtx )
        df_mtx = df_mtx.groupby( [ 'barcode', 'id_feature' ] ).sum( )
        df_mtx.reset_index( drop = False, inplace = True )

        # combine features
        df_feature = pd.concat( l_df_feature )
        df_feature.drop_duplicates( inplace = True )

        MTX_10X_Write( df_mtx, df_feature, dir_folder_output_mtx_10x )