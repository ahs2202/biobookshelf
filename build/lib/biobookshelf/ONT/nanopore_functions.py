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
        
def Gene_Read_Length( dir_file_bam, dir_file_gtf, return_list_of_read_length = False, thres_mapq = 30 ) :
    """ 
    # 2021-04-23 19:19:33 
    'dir_file_bam' : BAM file containing aligned nanopore reads
    'dir_file_gtf' : GTF file containing gene annotations of the reference genome to which the nanopore reads have been aligned
    'return_list_of_read_length' : return a dictionary containing list of read length for each name_gene
    'thres_mapq' : threshold for mapping quality
    
    return a dataframe summarizing aligned read length distribution for each gene   
    """
    # read gtf file and build interval tree
    df_gtf = GTF_Read( dir_file_gtf, parse_attr = True )
    df_gtf_gene = PD_Select( df_gtf, feature = 'gene' )
    dict_it = dict( )
    for gene_name, seqname, start, end in df_gtf_gene[ [ 'gene_name', 'seqname', 'start', 'end' ] ].values :
        if seqname not in dict_it :
            dict_it[ seqname ] = intervaltree.IntervalTree( )
        dict_it[ seqname ].addi( start, end, [ gene_name ] )

    # retrieve aligned length of nanopore read for each gene
    dict_length = dict( )
    with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
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

def Gene_10X_Adaptor( dir_file_bam, dir_file_gtf, thres_mapq = 60, float_error_rate = 0.15 ) :
    """ 
    # 2021-04-26 21:06:12 
    'dir_file_bam' : BAM file containing aligned nanopore reads
    'dir_file_gtf' : GTF file containing gene annotations of the reference genome to which the nanopore reads have been aligned
    'thres_mapq' : threshold for mapping quality
    'float_error_rate' : error rate for searching 10X adaptor sequences
    
    return a dataframe containing adaptor counts for each gene
    """
    # read gtf file and build interval tree
    df_gtf = GTF_Read( dir_file_gtf, parse_attr = True )
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
    with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
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
    return df