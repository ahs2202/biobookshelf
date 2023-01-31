# load internal module
from biobookshelf.main import *
import biobookshelf.PKG as PKG

def Guppy_Run_and_Combine_Output( path_folder_nanopore_sequencing_data = None, flag_barcoding_was_used = False, path_folder_output_fastq = None, id_flowcell = None, id_lib_prep = None, id_barcoding_kit = None, flag_use_cpu = True, int_n_threads = 18, flag_read_splitting = False ) :
    """
    # 2023-01-27 23:38:33 
    Run Guppy basecaller on the nanopore sequencing datafiles in the given folder 
    Automatically detect 'id_flowcell' and 'id_lib_prep' from the metadata saved in the folder (they can be manually set through arguments)
    
    'path_folder_nanopore_sequencing_data' : Run Guppy basecaller on the nanopore sequencing datafiles in the given folder. a list of nanopore sequencing datafolders of multiple runs of the same sample (due to failed runs, etc.) can be also given
    'flag_barcoding_was_used' : flag of whether a barcoding kit was used during sequencing
    'id_flowcell' : manually set 'id_flowcell' for guppy_bascaller run.
    'id_lib_prep' : manually set 'id_lib_prep' for guppy_bascaller run.
    'id_barcoding_kit' : manually set 'id_barcoding_kit' for guppy_bascaller run.
    'flag_use_cpu' : use CPU
    'int_n_threads' : number of CPU threads to use 
    """

    
    flag_entry_point = PKG.Detect_Entry_Point( 'biobook' ) # detect whether an entry point was used
    if flag_entry_point :
        parser = argparse.ArgumentParser( description = "Run Guppy basecaller on the nanopore sequencing datafiles in the given folder . This program has been developed by Hyunsu An (2021/06/03)." )
        parser.add_argument( 'path_folder_nanopore_sequencing_data', metavar = 'path_folder_nanopore_sequencing_data', type = str, nargs='+', help = '(Required) the nanopore sequencing data folder(s)' )
        parser.add_argument( "-b", "--flag_barcoding_was_used", help = "Set a flag indicating a barcoding kit was used during sequencing", action = 'store_true' )
        parser.add_argument( "-o", "--path_folder_output_fastq", help = "(optional) copy fastq files in the guppy output folder to the given directory" )
        parser.add_argument( "-F", "--id_flowcell", help = "(optional) explicitly define the flowcell type used in sequencing. e.g. FLO-MIN106" )
        parser.add_argument( "-L", "--id_lib_prep", help = "(optional) explicitly define the library sequencing kit used in sequencing. e.g. SQK-LSK109" )
        parser.add_argument( "-B", "--id_barcoding_kit", help = "(optional) explicitly define the library barcoding kit used in sequencing. e.g. EXP-NBD104" )
        parser.add_argument( "-C", "--flag_use_cpu", help = "(optional) use CPU only (use when GPU is not available)", action = 'store_true' )
        parser.add_argument( "-S", "--flag_read_splitting", help = "(optional) perform read splitting by using the Guppy 6's optional argument", action = 'store_true' )
        parser.add_argument( "-t", "--int_n_threads", help = "(default: 18) number of CPU threads to use", default = '18' )
        args = parser.parse_args( )

        # [input] parse arguments from parse_args
        path_folder_nanopore_sequencing_data = args.path_folder_nanopore_sequencing_data
        flag_barcoding_was_used = args.flag_barcoding_was_used
        path_folder_output_fastq = args.path_folder_output_fastq
        id_flowcell = args.id_flowcell
        id_lib_prep = args.id_lib_prep
        id_barcoding_kit = args.id_barcoding_kit
        flag_use_cpu = args.flag_use_cpu
        int_n_threads = int( args.int_n_threads )
        flag_read_splitting = args.flag_read_splitting
        
    
    ''' [parse arguments] '''
    if path_folder_nanopore_sequencing_data is None :
        print( "required arguments are not given, exiting" )
        if flag_entry_point :
            sys.exit( ) 
        else :
            return -1
    
    # process path_folder in 'l_path_folder_nanopore_sequencing_data'
    l_path_folder_nanopore_sequencing_data = path_folder_nanopore_sequencing_data if isinstance( path_folder_nanopore_sequencing_data, ( list ) ) else [ path_folder_nanopore_sequencing_data ] # set 'l_path_folder_nanopore_sequencing_data' according to the given 'path_folder_nanopore_sequencing_data'
    l = [ ]
    for path_folder in l_path_folder_nanopore_sequencing_data :
        path_folder = os.path.abspath( path_folder )
        if path_folder[ -1 ] != '/' : # add '/' to the end of the directory
            path_folder += '/'
        l.append( path_folder )
    l_path_folder_nanopore_sequencing_data = l
    # process 'path_folder_output_fastq'
    if path_folder_output_fastq is not None : # if output folder of fastq files was given
        path_folder_output_fastq = os.path.abspath( path_folder_output_fastq )
        if path_folder_output_fastq[ -1 ] != '/' : # add '/' to the end of the directory
            path_folder_output_fastq += '/'
            
            
    """ run Guppy basecaller for each nanopore sequencing data folder """
    set_path_file_fastq_gz = set( ) # a set of output fastq files
    for path_folder_nanopore_sequencing_data in l_path_folder_nanopore_sequencing_data : 
        # [input] parse arguments
        path_folder_nanopore_sequencing_data = os.path.abspath( path_folder_nanopore_sequencing_data )
        if path_folder_nanopore_sequencing_data[ -1 ] != '/' : # add '/' to the end of the directory
            path_folder_nanopore_sequencing_data += '/'

        # define and create a folder for all fast5 files
        path_folder_fast5 = path_folder_nanopore_sequencing_data + 'fast5/'
        if not os.path.exists( path_folder_fast5 ) : # if the folder already exists, skip moving fast5 files to a single folder (when fast-bascalling option has not been turned on, fast5 folder is present)
            os.makedirs( path_folder_fast5, exist_ok = True )
            # move all fast5 files into one folder
            l_path_file = glob.glob( f"{path_folder_nanopore_sequencing_data}*/*fast5" ) + glob.glob( f"{path_folder_nanopore_sequencing_data}*/*/*fast5" ) # retrieve fast5 files
            for path_file in l_path_file : # no barcoding
                shutil.copyfile( path_file, path_folder_fast5 + path_file.rsplit( '/', 1 )[ 1 ] )
        # retrieve flowcell type and library preperation method using summary file inside the directory

        if id_flowcell is None or id_lib_prep is None :
            l = glob.glob( f"{path_folder_nanopore_sequencing_data}final_summary_*.txt" )
            if len( l ) > 0 : 
                with open( l[ 0 ] ) as file :
                    content = file.read( ).strip( )
                    l_line = content.split( '\n' )
                if 'FLO-MIN114' in content : # 10.4.1 moter speed option
                    _, id_flowcell, id_lib_prep, str_moter_speed = list( line.split( 'protocol=' )[ 1 ].split( ':' ) for line in l_line if 'protocol=' == line[ : len( 'protocol=' ) ] )[ 0 ] # retrieve flowcell type and library preperation method using summary file inside the directory
                else :
                    _, id_flowcell, id_lib_prep = list( line.split( 'protocol=' )[ 1 ].split( ':' ) for line in l_line if 'protocol=' == line[ : len( 'protocol=' ) ] )[ 0 ] # retrieve flowcell type and library preperation method using summary file inside the directory
            
            if id_flowcell is None or id_lib_prep is None :
                l = glob.glob( f"{path_folder_nanopore_sequencing_data}report_*.md" ) # retry
                if len( l ) > 0 : 
                    with open( l[ 0 ] ) as file :
                        l_line = file.read( ).strip( ).split( '\n' )
                    _, id_flowcell, id_lib_prep = list( e.replace( '",', '' ) for e in list( line.split( '"exp_script_name":' )[ 1 ].split( ':' ) for line in l_line if '"exp_script_name":' in line )[ 0 ] ) # retrieve flowcell type and library preperation method using summary file inside the directory
                else :
                    print( "[Guppy_Run_and_Combine_Output] appropriate 'id_flowcell' and/or 'id_lib_prep' were not found, exiting" )
                    return -1 
        ''' run guppy basecaller and write output as a text file '''
        path_folder_guppy_output = f"{path_folder_nanopore_sequencing_data}guppy_out/"
        if not os.path.exists( path_folder_guppy_output ) :
            # compose guppy basecaller arguments
            l_args = [ 'guppy_basecaller' ] 
            if not flag_use_cpu : # use GPU if it exists
                l_args += [ '--device', 'auto' ]
            if flag_read_splitting :
                l_args += [ '--do_read_splitting' ] 
            l_args += [ '--cpu_threads_per_caller', str( int_n_threads ) ]

            if id_lib_prep == 'SQK-RBK096' : # change 'id_lib_prep' to guppy-compatible id
                id_lib_prep = 'SQK-RBK110-96'
            if flag_barcoding_was_used :
                ''' automatically set barcoding_kit '''
                if id_barcoding_kit is None :
                    id_barcoding_kit = "EXP-NBD104 EXP-NBD114" if 'LSK' in id_lib_prep else id_lib_prep
                print( ' '.join( l_args + [ "--flowcell", id_flowcell, "--kit", id_lib_prep, "--barcode_kits", id_barcoding_kit, "--compress_fastq", "--input_path", path_folder_fast5, "--save_path", path_folder_guppy_output ] ) )  
                run_guppy = subprocess.run( l_args + [ "--flowcell", id_flowcell, "--kit", id_lib_prep, "--barcode_kits", id_barcoding_kit, "--compress_fastq", "--input_path", path_folder_fast5, "--save_path", path_folder_guppy_output ], capture_output = True )
            else :
                run_guppy = subprocess.run( l_args + [ "--flowcell", id_flowcell, "--kit", id_lib_prep, "--compress_fastq", "--input_path", path_folder_fast5, "--save_path", path_folder_guppy_output ], capture_output = True )
            with open( path_folder_nanopore_sequencing_data + 'guppy_basecaller.out', 'w' ) as file :
                file.write( run_guppy.stdout.decode( ) )
        
        ''' combine fastq.gz output files of guppy_basecaller output '''
        if flag_barcoding_was_used :
            for path_folder_barcode in glob.glob( path_folder_guppy_output + 'pass/*/' ) + glob.glob( path_folder_guppy_output + 'fail/*/' ) :
                name_barcode = path_folder_barcode.rsplit( '/', 2 )[ 1 ] # retrieve barcode name from the path
                path_file_fastq_gz = f"{path_folder_guppy_output}{name_barcode}.fastq.gz"
                OS_FILE_Combine_Files_in_order( list( glob.glob( f"{path_folder_guppy_output}*/{name_barcode}/*" + '*fastq.gz' ) ), path_file_fastq_gz, flag_bgzip_output = False ) # combine fastq files
                set_path_file_fastq_gz.add( path_file_fastq_gz )
        else :
            path_file_fastq_gz = f"{path_folder_guppy_output}guppy_basecalled.fastq.gz"
            OS_FILE_Combine_Files_in_order( list( glob.glob( path_folder_guppy_output + '*fastq.gz' ) ), path_file_fastq_gz, flag_bgzip_output = False ) # combine fastq files
            set_path_file_fastq_gz.add( path_file_fastq_gz )
    ''' copy and combine output fastq files to the output directory '''
    if path_folder_output_fastq is not None : # if output folder of fastq files was given
        # group 'path_file_fastq_gz' based on 'name_file' 
        
        dict_name_file_to_dir = dict( )
        for d in set_path_file_fastq_gz :
            name_file = d.rsplit( '/', 1 )[ 1 ] 
            if name_file not in dict_name_file_to_dir :
                dict_name_file_to_dir[ name_file ] = [ ]
            dict_name_file_to_dir[ name_file ].append( d )
        
        for name_file in dict_name_file_to_dir : # for each output 'name_file' (barcodes), combine and copy the file to the given output folder
            l_path_file = dict_name_file_to_dir[ name_file ] # retrieve list of file path of the current barcode
            if len( l_path_file ) == 1 : # when only single file exist for the current barcode , simply copy the file to the output 
                OS_Run( [ 'cp', l_path_file[ 0 ], f"{path_folder_output_fastq}{name_file}" ] )
            else : # combine files 
                OS_FILE_Combine_Files_in_order( l_path_file, f"{path_folder_output_fastq}{name_file}", flag_bgzip_output = False ) # 
        
def Minimap2_Align( path_file_fastq, path_file_minimap2_index = '/node210data/shared/ensembl/Mus_musculus/index/minimap2/Mus_musculus.GRCm38.dna.primary_assembly.k_14.idx', path_folder_minimap2_output = None, n_threads = 20, verbose = True, drop_unaligned = False, return_bash_shellscript = False, n_threads_for_sort = 1, flag_use_split_prefix : bool = False ) :
    """ 
    # 2022-09-24 19:15:33 
    align given fastq file of nanopore reads using minimap2 and write an output as a bam file 
    'path_file_fastq' : input fastq or fasta file (gzipped or uncompressed file is accepted)
    'path_file_minimap2_index' : minimap2 index file
    'path_folder_minimap2_output' : minimap2 output folder
    'drop_unaligned' : a flag indicating whether reads not aligned to the reference ('SAM flag == 4') are included in the output bam file
    'return_bash_shellscript' : return shellscript instead of running minimap2 using the subprocess module
    'flag_use_split_prefix' = False # for large index, split-prefix should be used
    """
    path_folder_fastq, name_file_fastq = path_file_fastq.rsplit( '/', 1 )
    if path_folder_minimap2_output is None : # default output folder is a subdirectory of the folder containing the input fastq file
        path_folder_minimap2_output = f'{path_folder_fastq}/minimap2/'
    if path_folder_minimap2_output[ -1 ] != '/' : # add '/' at the end of the output directory if it does not exist
        path_folder_minimap2_output += '/'
    os.makedirs( path_folder_minimap2_output, exist_ok = True ) # create folder if it does not exist

    path_file_sam = f"{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.sam"
    path_file_bam = f"{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.bam"
    
    l_bash_shellscript = [ ]
    
    ''' perform minimap2 alignment '''
    l_arg = [ 'minimap2', '-t', str( int( n_threads ) ), '-ax', 'splice', "-o", path_file_sam ] 
    
    # for large index, split-prefix should be used
    if flag_use_split_prefix :
        l_arg += [ f'--split-prefix={path_folder_minimap2_output}{UUID( )}' ]
    
    if drop_unaligned :
        l_arg += [ '--sam-hit-only' ]
    l_arg += [ path_file_minimap2_index, path_file_fastq ]
    if return_bash_shellscript : # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
    else :
        run = subprocess.run( l_arg, capture_output = True )
        with open( f'{path_folder_minimap2_output}{name_file_fastq}.minimap2_aligned.out', 'w' ) as file :
            file.write( run.stdout.decode( ) )
        if verbose :
            print( 'minimap2 completed' )
            
    ''' sort output SAM file '''
    l_arg = [ 'samtools', 'sort', '-@', str( int( min( n_threads_for_sort, 10 ) ) ), '-O', "BAM", '-o', path_file_bam, path_file_sam ]
    if return_bash_shellscript : # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
        l_bash_shellscript.append( ' '.join( [ 'rm', '-f', path_file_sam ] ) )
    else :     
        run = subprocess.run( l_arg, capture_output = False )
        os.remove( path_file_sam ) # remove sam file
        
    ''' index resulting BAM file '''
    l_arg = [ 'samtools', 'index', path_file_bam ]
    if return_bash_shellscript : # perform minimap2 alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
    else :   
        run = subprocess.run( l_arg, capture_output = False )
        if verbose :
            print( 'samtools bam file compressing and indexing completed' )
    
    if return_bash_shellscript : # retrun bash shell scripts
        return ' && '.join( l_bash_shellscript )
    
def Minimap2_Index( path_file_fasta, path_file_minimap2_index = None, verbose = False ) :
    """ 
    # 2021-03-24 00:44:51 
    index given fasta file for nanopore reads alignment
    'path_file_fasta' : input reference fasta file
    'path_file_minimap2_index' : minimap2 index file
    """
    path_folder_fastq, name_file_fasta = path_file_fasta.rsplit( '/', 1 )
    if path_file_minimap2_index is None : # set the default directory of the minimap index 
        path_file_minimap2_index = f'{path_folder_fastq}/index/minimap2/{name_file_fasta}.ont.mmi'
    path_folder_minimap2_index, name_file_index = path_file_minimap2_index.rsplit( '/', 1 )
    path_folder_minimap2_index += '/'
    os.makedirs( path_folder_minimap2_index, exist_ok = True ) # create folder if it does not exist
    if os.path.exists( path_file_minimap2_index ) : # exit if an index file already exists
        return 
    # build minimap2 index
    run = subprocess.run( [ 'minimap2', '-x', 'map-ont', '-d', path_file_minimap2_index, path_file_fasta ], capture_output = True )
    
    with open( f'{path_folder_minimap2_index}{name_file_index}.minimap2_index.out', 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( 'minimap2 indexing completed' )
        
        
def Desalt_Index( path_file_fasta, path_folder_desalt_index = None, verbose = False ) :
    """ 
    # 2021-08-02 11:26:56 
    index given fasta file for nanopore reads alignment using deSALT
    'path_file_fasta' : input reference fasta file
    'path_folder_desalt_index' : a folder where the index will be saved.
    """
    path_folder_fastq, name_file_fasta = path_file_fasta.rsplit( '/', 1 )
    if path_folder_desalt_index is None : # set the default directory of the minimap index 
        os.makedirs( f'{path_folder_fastq}/index/desalt/', exist_ok = True )
        path_folder_desalt_index = f'{path_folder_fastq}/index/desalt/{name_file_fasta}.desalt.idx/'
    # build desalt index
    run = subprocess.run( [ 'deSALT', 'index', path_file_fasta, path_folder_desalt_index ], capture_output = True )
    with open( f"{path_folder_desalt_index.rsplit( '/', 1 )[ 0 ]}.desalt_index.out", 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( f'deSALT indexing of {path_file_fasta} completed' )
        
def Desalt_Align( path_file_fastq, path_folder_desalt_index = None, path_folder_desalt_output = None, arg_read_type = 'ont2d', path_file_gtf = None, n_threads = 20, verbose = True, return_bash_shellscript = False ) :
    """ 
    # 2021-06-14 23:19:19 
    align given fastq file of nanopore reads using deSALT and write an output as a bam file 
    'path_file_fastq' : input fastq or fasta file (gzipped or uncompressed file is accepted)
    'path_folder_desalt_index' : desalt index file
    'path_folder_desalt_output' : desalt output folder
    'arg_read_type' : parameter sets for each read type. available parameter sets are the following: 
                                   'ccs' (PacBio SMRT CCS reads): error rate 1%
                                   'clr' (PacBio SMRT CLR reads): error rate 15%
                                   'ont1d' (Oxford Nanopore 1D reads): error rate > 20%
                                   'ont2d' (Oxford Nanopore 2D reads): error rate > 12%
    'path_file_gtf' : directory of reference transcriptome annotations (for more accurate alignment)
    'return_bash_shellscript' : return shellscript instead of running desalt using the subprocess module
    """
    path_folder_fastq, name_file_fastq = path_file_fastq.rsplit( '/', 1 )
    if path_folder_desalt_output is None : # default output folder is a subdirectory of the folder containing the input fastq file
        path_folder_desalt_output = f'{path_folder_fastq}/desalt/'
    if path_folder_desalt_output[ -1 ] != '/' : # add '/' at the end of the output directory if it does not exist
        path_folder_desalt_output += '/'
    os.makedirs( path_folder_desalt_output, exist_ok = True ) # create folder if it does not exist

    path_file_sam = f"{path_folder_desalt_output}{name_file_fastq}.desalt_aligned.sam"
    path_file_bam = f"{path_folder_desalt_output}{name_file_fastq}.desalt_aligned.bam"
    
    l_bash_shellscript = [ ]
    
    ''' perform desalt alignment '''
    l_arg = [ 'deSALT', 'aln', '-t', str( int( n_threads ) ), '-x', arg_read_type, "-o", path_file_sam ] 
    if path_file_gtf is not None :
        path_file_gtf_info = f"{path_file_gtf}.desalt.info" # directory of a processed gtf file
        # process a given gtf file into deSALT GTF info
        if not os.path.exists( path_file_gtf_info ) :
            run = subprocess.run( [ 'Annotation_Load.py', path_file_gtf, path_file_gtf_info ], capture_output = True )
            with open( f'{path_file_gtf}.desalt_annotation_load.out', 'w' ) as file :
                file.write( run.stdout.decode( ) )
        l_arg += [ '-G', path_file_gtf ]
    l_arg += [ path_folder_desalt_index, path_file_fastq ]
    if return_bash_shellscript : # perform desalt alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
    else :
        run = subprocess.run( l_arg, capture_output = True )
        with open( f'{path_folder_desalt_output}{name_file_fastq}.desalt_aligned.out', 'w' ) as file :
            file.write( run.stdout.decode( ) )
        if verbose :
            print( 'desalt completed' )
            
    ''' sort output SAM file '''
    l_arg = [ 'samtools', 'sort', '-@', str( int( min( n_threads, 10 ) ) ), '-O', "BAM", '-o', path_file_bam, path_file_sam ]
    if return_bash_shellscript : # perform desalt alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
        l_bash_shellscript.append( ' '.join( [ 'rm', '-f', path_file_sam ] ) )
    else :     
        run = subprocess.run( l_arg, capture_output = False )
        os.remove( path_file_sam ) # remove sam file
        
    ''' index resulting BAM file '''
    l_arg = [ 'samtools', 'index', path_file_bam ]
    if return_bash_shellscript : # perform desalt alignment using subprocess module
        l_bash_shellscript.append( ' '.join( l_arg ) )
    else :   
        run = subprocess.run( l_arg, capture_output = False )
        if verbose :
            print( 'samtools bam file compressing and indexing completed' )
    
    if return_bash_shellscript : # retrun bash shell scripts
        return ' && '.join( l_bash_shellscript )
        
def FeatureCounts( path_file_annotation, path_file_output, * l_path_file_input, str_type_attribute = 'gene_id', verbose = False, return_dataframe = True ) :
    """ 
    # 2021-03-07 18:44:04 
    count aligned reads using featureCounts
    'str_type_attribute' : see featureCounts help messages
    'return_dataframe' : return a dataframe containing the featureCounts output
    """
    path_folder_output, name_file_output = path_file_output.rsplit( '/', 1 )
    path_folder_output += '/'
    os.makedirs( path_folder_output, exist_ok = True ) # create folder if it does not exist
    
    run = subprocess.run( [ 'featureCounts', '-L', '-a', path_file_annotation, '-o', path_file_output, '-g', str_type_attribute ] + list( l_path_file_input ), capture_output = True )
    with open( f'{path_file_output}.featureCounts.out', 'w' ) as file :
        file.write( run.stdout.decode( ) )
    if verbose :
        print( 'featureCounts completed' )
    if return_dataframe : # return a dataframe containing the featureCounts output
        df = pd.read_csv( path_file_output, sep = '\t', low_memory = False, skiprows = [ 0 ] )
        return df

def Gene_Read_Count( path_file_bam, dict_it_gene, thres_mapq = 60 ) :
    """ 
    Count number of reads uniquely aligned to each gene 
    
    'path_file_bam' : BAM file containing aligned nanopore reads
    'dict_it_gene' : dictionary of interval trees from a GTF file containing gene annotations of the reference genome to which the nanopore reads have been aligned (returned by 'GTF_Interval_Tree')
    'thres_mapq' : threshold for mapping quality
    
    return a dataframe
    """

    # retrieve aligned length of nanopore read for each gene
    dict_count = dict( )
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            if r.reference_name not in dict_it_gene : # skip read aligned to segment that does not contain genes
                continue
            if r.mapq < thres_mapq : # skip read whose mapq is below 'thres_mapq'
                continue
            if r.seq is None : # skip multi-mapped reads (minimap2 skip sequence information for multi-mapped reads)
                continue
            set_overlap = dict_it_gene[ r.reference_name ].overlap( r.reference_start, r.reference_end )
            if len( set_overlap ) == 0 : # skip if the read has not been aligned to any genes
                continue
            name_gene = list( set_overlap )[ 0 ][ 2 ][ 0 ]
            if name_gene not in dict_count :
                dict_count[ name_gene ] = 0
            dict_count[ name_gene ] += 1
    return dict_count
    
    # summarize read length distribution for each gene
    l_l = [ ]
    for name_gene in dict_length :
        arr = np.array( dict_length[ name_gene ], dtype = int )
        l_l.append( [ name_gene, len( arr ), np.mean( arr ), np.std( arr ) ] )
    df_gene_read_length_summary = pd.DataFrame( l_l, columns = [ 'name_gene', 'n_reads', 'mean_read_length', 'std_read_length' ] )
    return df_gene_read_length_summary
    
def Gene_Read_Length( path_file_bam, path_file_gtf, return_list_of_read_length = False, thres_mapq = 30 ) :
    """ 
    # 2021-04-23 19:19:33 
    'path_file_bam' : BAM file containing aligned nanopore reads
    'path_file_gtf' : GTF file containing gene annotations of the reference genome to which the nanopore reads have been aligned
    'return_list_of_read_length' : return a dictionary containing list of read length for each name_gene
    'thres_mapq' : threshold for mapping quality
    
    return a dataframe summarizing aligned read length distribution for each gene   
    """
    # read gtf file and build interval tree
    df_gtf = GTF_Read( path_file_gtf, parse_attr = True )
    df_gtf_gene = PD_Select( df_gtf, feature = 'gene' )
    dict_it = dict( )
    for gene_name, seqname, start, end in df_gtf_gene[ [ 'gene_name', 'seqname', 'start', 'end' ] ].values :
        if seqname not in dict_it :
            dict_it[ seqname ] = intervaltree.IntervalTree( )
        dict_it[ seqname ].addi( start, end, [ gene_name ] )

    # retrieve aligned length of nanopore read for each gene
    dict_length = dict( )
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            if r.reference_name not in dict_it : # skip read aligned to segment that does not contain genes
                continue
            if r.mapq < thres_mapq : # skip read whose mapq is below 'thres_mapq'
                continue
            set_overlap = dict_it[ r.reference_name ].overlap( r.reference_start, r.reference_end )
            if len( set_overlap ) == 0 :
                continue
            name_gene = list( set_overlap )[ 0 ][ 2 ][ 0 ]
            if name_gene not in dict_length :
                dict_length[ name_gene ] = [ ]
            dict_length[ name_gene ].append( r.qlen )
    if return_list_of_read_length : # if 'return_list_of_read_length' is set to True, return a dictionary containing list of read length for each name_gene 
        return dict_length
    
    # summarize read length distribution for each gene
    l_l = [ ]
    for name_gene in dict_length :
        arr = np.array( dict_length[ name_gene ], dtype = int )
        l_l.append( [ name_gene, len( arr ), np.mean( arr ), np.std( arr ) ] )
    df_gene_read_length_summary = pd.DataFrame( l_l, columns = [ 'name_gene', 'n_reads', 'mean_read_length', 'std_read_length' ] )
    return df_gene_read_length_summary

def Gene_10X_Adaptor( path_file_bam, path_file_gtf, thres_mapq = 60, float_error_rate = 0.15 ) :
    """ 
    # 2021-04-27 11:30:13 
    'path_file_bam' : BAM file containing aligned nanopore reads
    'path_file_gtf' : GTF file containing gene annotations of the reference genome to which the nanopore reads have been aligned
    'thres_mapq' : threshold for mapping quality
    'float_error_rate' : error rate for searching 10X adaptor sequences
    
    return a dataframe containing adaptor counts for each gene
    """
    # read gtf file and build interval tree
    df_gtf = GTF_Read( path_file_gtf, parse_attr = True )
    df_gtf_gene = PD_Select( df_gtf, feature = 'gene' )
    dict_it = dict( )
    for gene_name, seqname, start, end in df_gtf_gene[ [ 'gene_name', 'seqname', 'start', 'end' ] ].values :
        if seqname not in dict_it :
            dict_it[ seqname ] = intervaltree.IntervalTree( )
        dict_it[ seqname ].addi( start, end, [ gene_name ] )

    # adaptor sequences for counting adaptor sequences from the soft-clipped sequences at the ends
    str_seq_tso = "AAGCAGTGGTATCAACGCAGAGTACAT"
    str_seq_r1 = "CTACACGACGCTCTTCCGATCT"
    str_seq_tso_rc = NGS_SEQ_Reverse_Complement( "AAGCAGTGGTATCAACGCAGAGTACAT" )
    str_seq_r1_rc = NGS_SEQ_Reverse_Complement( "CTACACGACGCTCTTCCGATCT" )

    dict_adaptor = dict( )
    with pysam.AlignmentFile( path_file_bam, 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            if r.reference_name not in dict_it : # skip read aligned to segment that does not contain genes
                continue
            if r.mapq < thres_mapq : # skip read whose mapq is below 'thres_mapq'
                continue
            if r.seq is None : # skip multi-mapped reads (minimap2 skip sequence information for multi-mapped reads)
                continue
            set_overlap = dict_it[ r.reference_name ].overlap( r.reference_start, r.reference_end )
            if len( set_overlap ) == 0 : # skip if no gene was assigned to the read 
                continue
            name_gene = list( set_overlap )[ 0 ][ 2 ][ 0 ] # retrieve name of the gene that the read was assigned

            # retrieve soft-clipped sequences from the ends of the read
            l_cigar = list( NGS_CIGAR_Iterate_a_CIGAR_String( r.cigarstring ) ) # retrieve tuples of a cigar string
            l_seq = [ ] 
            if l_cigar[ 0 ][ 1 ] == 'S' :
                l_seq.append( r.seq[ : l_cigar[ 0 ][ 0 ] ] )
            if l_cigar[ -1 ][ 1 ] == 'S' :
                l_seq.append( r.seq[ - l_cigar[ -1 ][ 0 ] : ] )
            # count adaptors at the soft-clipped sequences at the ends
            n_tso_count, n_r1_count = 0, 0
            for seq_adaptor in [ str_seq_tso, str_seq_tso_rc ] :
                for seq in l_seq :
                    if STR.Search_Subsequence( seq, seq_adaptor, error_rate = float_error_rate )[ 'matched_subsequence' ] is not None : # increase the count if the adaptor is found in the soft-clipped bsequence
                        n_tso_count += 1
            for seq_adaptor in [ str_seq_r1, str_seq_r1_rc ] :
                for seq in l_seq :
                    if STR.Search_Subsequence( seq, seq_adaptor, error_rate = float_error_rate )[ 'matched_subsequence' ] is not None :
                        n_r1_count += 1

            # count each adaptor_type (classification of reads based on adaptor counts) for each gene
            if name_gene not in dict_adaptor :
                dict_adaptor[ name_gene ] = dict( )
            name_adaptor_type = f'tso_{n_tso_count}__r1_{n_r1_count}'
            if name_adaptor_type not in dict_adaptor[ name_gene ] :
                dict_adaptor[ name_gene ][ name_adaptor_type ] = 0
            dict_adaptor[ name_gene ][ name_adaptor_type ] += 1
    
    # compose a dataframe
    df = pd.DataFrame( dict_adaptor ).T
    df.fillna( 0, inplace = True )
    df[ 'total_count' ] = df.sum( axis = 1 )
    df.sort_values( 'total_count', ascending = False, inplace = True )
    df.index.name = 'name_gene'
    return df


def Check_plasmid_with_nanopore_sequencing( path_file_fasta_ref = None, path_file_fastq = None, path_folder_output = 'default', n_threads = 10, flag_correct_reference = False, flag_perform_assembly = False ) :
    """
    # 2021-06-03 21:52:00 
    Correct plasmid sequence using nanopore sequencing 
    
    'path_file_fastq' : (Required) directory of a fasta file containing the reference sequence. (e.g. a fasta sequence from AddGene)
    'path_folder_output' : (Default: subdirectory of the folder containing the given fastq file) directory of output folder
    'n_threads' : number of threads to run the pipeline
    
    """
    # 'flag_correct_reference' : manually set 'id_barcoding_kit' for guppy_bascaller run.
    
    flag_entry_point = PKG.Detect_Entry_Point( 'biobook' ) # detect whether an entry point was used
    if flag_entry_point :
        parser = argparse.ArgumentParser( description = "Analyze nanopore sequencing data of a plasmid. This program has been developed by Hyunsu An (2021/06/03)." )
        parser.add_argument( "-r", "--path_file_fasta_ref", help = "(Required) directory of a fasta file containing the reference sequence. (e.g. a fasta sequence from AddGene)" )
        parser.add_argument( "-i", "--path_file_fastq", help = "(Required) directory of a fastq file from a nanopore sequencing" )
        parser.add_argument( "-o", "--path_folder_output", help = "(Default: subdirectory of the folder containing the given fastq file) directory of output folder", default = 'default' )
        parser.add_argument( "-t", "--threads", help = "(Default: 10) Number of threads to use in the current compute node.", default = '10' )
        parser.add_argument( "-a", "--flag_perform_assembly", help = "(Default: False) Perform assembly using Flye", action = 'store_true' )

        args = parser.parse_args( )
        # [input] parse arguments from parse_args
        path_file_fasta_ref = args.path_file_fasta_ref
        path_file_fastq = args.path_file_fastq
        path_folder_output = args.path_folder_output
        flag_perform_assembly = args.flag_perform_assembly
        # flag_correct_reference = args.flag_correct_reference
        n_threads = int( args.threads )
    
    
    ''' [parse arguments] '''
    if path_file_fastq is None or path_file_fastq is None  :
        print( "required arguments are not given, exiting" )
        if flag_entry_point :
            sys.exit( ) 
        else :
            return -1
        

    # [process input] output folder 
    # retrieve absolute paths
    path_file_fasta_ref = os.path.abspath( path_file_fasta_ref )
    path_file_fastq = os.path.abspath( path_file_fastq )
    # set default output folder 
    if path_folder_output == 'default' : 
        path_folder, name_file = path_file_fastq.rsplit( '/', 1 )
        path_folder += '/'
        if name_file.rsplit( '.', 1 )[ 1 ].lower( ) == 'gz' :
            name_file = name_file.rsplit( '.', 1 )[ 0 ]
        name_file = name_file.rsplit( '.', 1 )[ 0 ]
        path_folder_output = f"{path_folder}{name_file}/"
    path_folder_output = os.path.abspath( path_folder_output )
    if path_folder_output[ -1 ] != '/' : # add '/' at the end of the output folder
        path_folder_output += '/'
    if os.path.exists( path_folder_output ) : # if output folder already exists, exit
        print( 'output folder already exists, exiting' )
        if flag_entry_point :
            sys.exit( ) 
        else :
            return -1
    else : # create an output folder
        os.makedirs( path_folder_output )
            
    """ analyze plasmid sequence """
    path_file_index_ref = f"{path_folder_output}index.ont.mmi"
    Minimap2_Index( path_file_fasta_ref, path_file_index_ref )
    path_folder_minimap2_output = f"{path_folder_output}minimap2/"
    Minimap2_Align( path_file_fastq, path_file_minimap2_index = path_file_index_ref, path_folder_minimap2_output = path_folder_minimap2_output, n_threads = n_threads )

    dict_fasta = FASTA_Read( path_file_fasta_ref )
    str_fasta_ref = dict_fasta[ list( dict_fasta )[ 0 ] ].upper( ) # read genome sequence of the reference (first sequence of the reference fasta file)
    shutil.copyfile( path_file_fasta_ref, f"{path_folder_minimap2_output}{path_file_fasta_ref.rsplit( '/', 1 )[ 1 ]}" ) # copy reference sequence file to minimap2 output file (convenient minimap visualization

    # count number of each base for each position of the reference genome
    str_list_of_bases = "ATGC-"
    dict_base_to_index = dict( ( ( str_base, i ) for i, str_base in enumerate( str_list_of_bases ) ) )
    arr_count_base = np.zeros( ( 5, len( str_fasta_ref ) ) )

    l = [ ]
    with pysam.AlignmentFile( glob.glob( f"{path_folder_minimap2_output}*.bam" )[ 0 ], 'rb' ) as samfile :
        for r in samfile.fetch( ) :
            for pos_read, pos_ref in r.aligned_pairs :
                if pos_ref is not None and r.seq is not None : # ignore insertion mutation type # ignore secondary alignment (where sequence is None)
                    arr_count_base[ dict_base_to_index[ '-' if pos_read is None else r.seq[ pos_read ] ] ][ pos_ref ] += 1 # detect deletion by checking aligned positions of a read

    # summerize results
    arr_coverage = arr_count_base.sum( axis = 0 )
    arr_coverage[ arr_coverage == 0 ] = -1 # mask positions with zero coverage
    df_summary = pd.DataFrame( arr_count_base / arr_coverage, index = list( str_list_of_bases ), columns = np.arange( 1, len( str_fasta_ref ) + 1 ) ) # 1-based coordinates
    df_summary.loc[ 'coverage' ] = arr_count_base.sum( axis = 0 )
    str_fasta_consensus = ''.join( list( str_list_of_bases[ arr_freq.argmax( ) ] if arr_freq.sum( ) > 0.1 else '*' for arr_freq in ( arr_count_base / arr_coverage ).T ) ) # if arr_frequency contains only zero values (no coverage), put '*' in the sequence
    df_summary.loc[ 'consensus_sequence' ] = list( str_fasta_consensus )
    df_summary.loc[ 'reference_sequence' ] = list( str_fasta_ref )
    l = [ ]
    # classify mutations
    for base_consensus, base_ref in zip( df_summary.loc[ 'consensus_sequence' ].values, df_summary.loc[ 'reference_sequence' ].values ) :
        str_mut_type = '.' # default mut_type
        if base_consensus == base_ref :
            pass
        elif base_consensus == '-' :
            str_mut_type = 'deletion'
        elif base_ref == 'N' :
            pass
        elif base_consensus == '*' :
            str_mut_type = 'unknown'
        else :
            str_mut_type = 'substitution'
        l.append( str_mut_type )
    df_summary.loc[ 'mutation_type' ] = l
    df_summary.to_excel( f"{path_folder_output}summary__base_frequency.xlsx" )

    # detect substitutions and extract a flanking sequence for each substitution
    df_substitution = df_summary.T[ df_summary.T.mutation_type == 'substitution' ]
    df_substitution.index.name = 'location_of_mutation'
    df_substitution.reset_index( inplace = True )

    len_flanking_base = 5 # number of bases flanking a mutation
    l_l = list( )
    for index in df_substitution.location_of_mutation.values - 1 : # 0-based coordinates
        sl = slice( index - len_flanking_base, index + len_flanking_base + 1 )
        l_l.append( [ str_fasta_consensus[ sl ], str_fasta_ref[ sl ] ] )
    df_substitution = df_substitution.join( pd.DataFrame( l_l, columns = [ 'flanking_sequence_consensus', 'flanking_sequence_reference' ] ) )
    df_substitution.to_excel( f"{path_folder_output}summary__substitution.xlsx" )
    
    if flag_perform_assembly : # perform assembly if the flag is on
        pass
        
    
    
    
    