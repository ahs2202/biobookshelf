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
    
def MTX_10X_Combine( dir_folder_output_mtx_10x, * l_dir_folder_input_mtx_10x ) :
    '''
    # 2021-11-24 13:13:50 
    Combine 10X count matrix files from the given list of folders and write combined output files to the given output folder 'dir_folder_output_mtx_10x'
    '''
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