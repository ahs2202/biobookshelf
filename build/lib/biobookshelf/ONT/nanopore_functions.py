# load internal module
from biobookshelf.main import *

def Guppy_Run_and_Combine_Output( dir_folder_nanopore_sequencing_data, flag_barcoding_was_used = False, dir_folder_output_fastq = None ) :
    """
    # 2021-03-25 22:52:55 
    Run Guppy basecaller on the nanopore sequencing datafiles in the given folder 
    'dir_folder_nanopore_sequencing_data' : Run Guppy basecaller on the nanopore sequencing datafiles in the given folder 
    'flag_barcoding_was_used' : flag of whether a barcoding kit was used during sequencing
    """
    # [input] parse arguments
    dir_folder_nanopore_sequencing_data = Program__Get_Absolute_Path_of_a_File( dir_folder_nanopore_sequencing_data )
    if dir_folder_nanopore_sequencing_data[ -1 ] != '/' : # add '/' to the end of the directory
        dir_folder_nanopore_sequencing_data += '/'
    
    # define and create a folder for all fast5 files
    dir_folder_fast5 = dir_folder_nanopore_sequencing_data + 'fast5/'
    if os.path.exists( dir_folder_fast5 ) : # if the folder already exists, empty the folder
        shutil.rmtree( dir_folder_fast5 )
    os.makedirs( dir_folder_fast5, exist_ok = True )
    # move all fast5 files into one folder
    l_dir_file = glob.glob( f"{dir_folder_nanopore_sequencing_data}*/*fast5" ) + glob.glob( f"{dir_folder_nanopore_sequencing_data}*/*/*fast5" ) # retrieve fast5 files
    for dir_file in l_dir_file : # no barcoding
        shutil.copyfile( dir_file, dir_folder_fast5 + dir_file.rsplit( '/', 1 )[ 1 ] )
    # retrieve flowcell type and library preperation method using summary file inside the directory
    
    l = glob.glob( f"{dir_folder_nanopore_sequencing_data}final_summary_*.txt" )
    if len( l ) == 0 : 
        print( '[Guppy_Run_and_Combine_Output] no summary file found, exiting' )
        return -1
    with open( l[ 0 ] ) as file :
        l_line = file.read( ).strip( ).split( '\n' )
    _, id_flowcell, id_lib_prep = list( line.split( 'protocol=' )[ 1 ].split( ':' ) for line in l_line if 'protocol=' == line[ : len( 'protocol=' ) ] )[ 0 ] # retrieve flowcell type and library preperation method using summary file inside the directory
    
    # run guppy basecaller and write output as a text file
    dir_folder_guppy_output = f"{dir_folder_nanopore_sequencing_data}guppy_out/"
    if flag_barcoding_was_used :
        run_guppy = subprocess.run( [ 'guppy_basecaller', '--device', 'auto', '--cpu_threads_per_caller', '18', "--flowcell", id_flowcell, "--kit", id_lib_prep, "--barcode_kits", "EXP-NBD104" if 'LSK' in id_lib_prep else "SQK-RBK004", "--compress_fastq", "--input_path", dir_folder_fast5, "--save_path", dir_folder_guppy_output ], capture_output = True )
    else :
        run_guppy = subprocess.run( [ 'guppy_basecaller', '--device', 'auto', '--cpu_threads_per_caller', '18', "--flowcell", id_flowcell, "--kit", id_lib_prep, "--compress_fastq", "--input_path", dir_folder_fast5, "--save_path", dir_folder_guppy_output ], capture_output = True )
    with open( dir_folder_nanopore_sequencing_data + 'guppy_basecaller.out', 'w' ) as file :
        file.write( run_guppy.stdout.decode( ) )
    # combine fastq.gz output files of guppy_basecaller output
    l_dir_file_fastq_gz = [ ] # list of output fastq files
    if flag_barcoding_was_used :
        for dir_folder_barcode in glob.glob( dir_folder_guppy_output + '*/' ) :
            name_barcode = dir_folder_barcode.rsplit( '/', 2 )[ 1 ] # retrieve barcode name from the path
            dir_file_fastq_gz = f"{dir_folder_guppy_output}{name_barcode}.fastq.gz"
            OS_FILE_Combine_Files_in_order( glob.glob( dir_folder_barcode + '*fastq.gz' ), dir_file_fastq_gz, overwrite_existing_file = True )
            l_dir_file_fastq_gz.append( dir_file_fastq_gz )
    else :
        dir_file_fastq_gz = f"{dir_folder_guppy_output}guppy_basecalled.fastq.gz"
        OS_FILE_Combine_Files_in_order( glob.glob( dir_folder_guppy_output + '*fastq.gz' ), dir_file_fastq_gz, overwrite_existing_file = True )
        l_dir_file_fastq_gz.append( dir_file_fastq_gz )
        
    if dir_folder_output_fastq is not None : # if output folder of fastq files was given
        dir_folder_output_fastq = Program__Get_Absolute_Path_of_a_File( dir_folder_output_fastq )
        if dir_folder_output_fastq[ -1 ] != '/' : # add '/' to the end of the directory
            dir_folder_output_fastq += '/'
        for dir_file_fastq_gz in l_dir_file_fastq_gz : # for each output fastq file, copy the file to the given fastq folder
            shutil.copyfile( dir_file_fastq_gz, f"{dir_folder_output_fastq}{dir_file_fastq_gz.rsplit( '/', 1 )[ 1 ]}" )
            
            
        
def Minimap2_Align( dir_file_fastq, dir_file_minimap2_index = '/node210data/shared/ensembl/Mus_musculus/index/minimap2/Mus_musculus.GRCm38.dna.primary_assembly.k_14.idx', dir_folder_minimap2_output = None, n_threads = 20, verbose = True ) :
    """ 
    # 2021-03-06 21:37:21 
    align given fastq file of nanopore reads using minimap2 and write an output as a bam file 
    'dir_file_fastq' : input fastq or fasta file (gzipped or uncompressed file is accepted)
    'dir_file_minimap2_index' : minimap2 index file
    'dir_folder_minimap2_output' : minimap2 output folder
    """
    dir_folder_fastq, name_file_fastq = dir_file_fastq.rsplit( '/', 1 )
    if dir_folder_minimap2_output is None : # default output folder is a subdirectory of the folder containing the input fastq file
        dir_folder_minimap2_output = f'{dir_folder_fastq}/minimap2/'
    if dir_folder_minimap2_output[ -1 ] != '/' : # add '/' at the end of the output directory if it does not exist
        dir_folder_minimap2_output += '/'
    os.makedirs( dir_folder_minimap2_output, exist_ok = True ) # create folder if it does not exist

    dir_file_sam = f"{dir_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.sam"
    dir_file_bam = f"{dir_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.bam"
    # perform minimap2 alignment
    run = subprocess.run( [ 'minimap2', '-t', str( int( n_threads ) ), '-ax', 'splice', "-o", dir_file_sam, dir_file_minimap2_index, dir_file_fastq ], capture_output = True )
    with open( f'{dir_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.out', 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( 'minimap2 completed' )
    run = subprocess.run( [ 'samtools', 'sort', '-@', str( int( min( n_threads, 10 ) ) ), '-O', "BAM", '-o', dir_file_bam, dir_file_sam ], capture_output = False )
    run = subprocess.run( [ 'samtools', 'index', dir_file_bam ], capture_output = False )
    if verbose :
        print( 'samtools bam file compressing and indexing completed' )
        
def Minimap2_Index( dir_file_fasta, dir_file_minimap2_index = None, verbose = False ) :
    """ 
    # 2021-03-24 00:44:51 
    index given fasta file for nanopore reads alignment
    'dir_file_fasta' : input reference fasta file
    'dir_file_minimap2_index' : minimap2 index file
    """
    dir_folder_fastq, name_file_fasta = dir_file_fasta.rsplit( '/', 1 )
    if dir_file_minimap2_index is None : # set the default directory of the minimap index 
        dir_file_minimap2_index = f'{dir_folder_fastq}/index/minimap2/{name_file_fasta}.ont.mmi'
    dir_folder_minimap2_index, name_file_index = dir_file_minimap2_index.rsplit( '/', 1 )
    dir_folder_minimap2_index += '/'
    os.makedirs( dir_folder_minimap2_index, exist_ok = True ) # create folder if it does not exist
    if os.path.exists( dir_file_minimap2_index ) : # exit if an index file already exists
        return 
    # build minimap2 index
    run = subprocess.run( [ 'minimap2', '-x', 'map-ont', '-d', dir_file_minimap2_index, dir_file_fasta ], capture_output = True )
    
    with open( f'{dir_folder_minimap2_index}{name_file_index}.minimap2_index.out', 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( 'minimap2 indexing completed' )
        
def FeatureCounts( dir_file_annotation, dir_file_output, * l_dir_file_input, str_type_attribute = 'gene_id', verbose = False, return_dataframe = True ) :
    """ 
    # 2021-03-07 18:44:04 
    count aligned reads using featureCounts
    'str_type_attribute' : see featureCounts help messages
    'return_dataframe' : return a dataframe containing the featureCounts output
    """
    dir_folder_output, name_file_output = dir_file_output.rsplit( '/', 1 )
    dir_folder_output += '/'
    os.makedirs( dir_folder_output, exist_ok = True ) # create folder if it does not exist
    
    run = subprocess.run( [ 'featureCounts', '-L', '-a', dir_file_annotation, '-o', dir_file_output, '-g', str_type_attribute ] + list( l_dir_file_input ), capture_output = True )
    with open( f'{dir_file_output}.featureCounts.out', 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( 'featureCounts completed' )
    if return_dataframe : # return a dataframe containing the featureCounts output
        df = pd.read_csv( dir_file_output, sep = '\t', low_memory = False, skiprows = [ 0 ] )
        return df
        
