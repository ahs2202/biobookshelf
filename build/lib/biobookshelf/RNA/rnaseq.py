from biobookshelf.main import *

def STAR_Index( dir_file_fasta, dir_file_annotation, dir_folder_index = None, int_sjdboverhang = 75, n_threads = 10 ) :
    """
    # 2021-03-23 22:08:37 hyunsu-an
    create index of STAR aligner
    
    'dir_file_fasta' : genome fasta file
    'dir_file_annotation' : GTF file containing gene annotation
    'int_sjdboverhang' : see documentation of STAR
    """
    # set default 'dir_folder_index' 
    dir_folder_fasta, name_file_fasta = dir_file_fasta.rsplit( '/', 1 )
    if dir_folder_index is None : # set default output folder
        dir_folder_index = f'{dir_folder_fasta}/index/STAR/{name_file_fasta}_sjdbOverhang_{int( int_sjdboverhang )}/'  
    if dir_folder_index[ -1 ] != '/' : # add '/' at the end of the output directory if it does not exist
        dir_folder_index += '/'
    os.makedirs( dir_folder_index, exist_ok = True ) # create folder if it does not exist
    
    if os.path.exists( f"{dir_folder_index}SA" ) : # exit if STAR index exists
        print( f"[STAR_Index] it appears that STAR index at {dir_folder_index} already exists" )
        return -1
    
    run_star = subprocess.run( [ 'STAR', '--runThreadN', str( int( n_threads ) ), '--runMode', 'genomeGenerate', '--genomeDir', dir_folder_index, '--genomeFastaFiles', dir_file_fasta, '--sjdbGTFfile', dir_file_annotation, '--sjdbOverhang', str( int( int_sjdboverhang ) ) ], capture_output = True )
    
    with open( f'{dir_folder_index}{name_file_fasta}.STAR_index.out', 'w' ) as file :
        file.write( run_star.stdout.decode( ) )
        
def STAR_Align( dir_folder_index, dir_file_read_1, dir_file_read_2 = None, dir_prefix_output = None, n_threads = 10, index_bam = True ) :
    """
    # 2021-03-23 23:05:37 
    Align reads with STAR aligner
    
    'dir_folder_index' : dir_folder containing STAR index
    'dir_file_read_1' : file containing read_1 (gzipped file supported)
    'dir_file_read_2' : file containing read_2 (gzipped file supported). read_1 file and read_2 file should be all gzipped or unzipped.
    'dir_prefix_output' : prefix of the STAR output files
    'index_bam' : index bam file using samtools
    
    return 'dir_file_bam' : aligned BAM file 
    """
    # set n_threads
    n_threads = min( n_threads, 16 ) # max 16 threads
    # set default 'dir_prefix_output' 
    dir_folder_read, name_file_read = dir_file_read_1.rsplit( '/', 1 )
    dir_folder_read += '/'
    if dir_prefix_output is None : # set default prefix for STAR output
        os.makedirs( f'{dir_folder_read}STAR_out/', exist_ok = True ) # create an output folder if it does not exist
        dir_prefix_output = f'{dir_folder_read}STAR_out/{name_file_read}'  
    
    l_args = [ 'STAR', '--runMode', "alignReads", '--genomeDir', dir_folder_index, '--runThreadN', str( int( n_threads ) ), '--outFileNamePrefix', dir_prefix_output, '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMunmapped', 'Within', '--outSAMattributes', 'Standard' ] 
    if dir_file_read_1.rsplit( '.', 1 )[ 1 ].lower( ) == 'gz' : # if given input file containing reads is gzipped
        l_args += [ "--readFilesCommand", "zcat" ]
    # add read_1 and read_2 file to arguments
    l_args += [ "--readFilesIn", dir_file_read_1 ]
    if dir_file_read_2 is not None :
        l_args += [ dir_file_read_2 ]
    run_star = subprocess.run( l_args, capture_output = True )
    
    with open( f'{dir_prefix_output}.STAR_Align.out', 'w' ) as file :
        file.write( run_star.stdout.decode( ) )
    dir_file_bam = f"{dir_prefix_output}Aligned.sortedByCoord.out.bam"
    if index_bam : # index output bam file if 'index_bam' is True
        pysam.index( dir_file_bam )
    return dir_file_bam

def Coverage( dir_file_bam, thres_mapq = 200 ) :
    """
    # 2021-03-24 00:57:51 
    calculate coverage from BAM file 
    
    'thres_mapq' : threshold for MAPQ for calculating coverage
    return 'dict_refname_to_cov'
    """
    with pysam.AlignmentFile( dir_file_bam, 'rb' ) as samfile :
        sam_header = samfile.header
        dict_refname_to_cov = dict( ( r[ 'SN' ], np.zeros( r[ 'LN' ], dtype = int ) ) for r in sam_header.to_dict( )[ 'SQ' ] ) # dictionary containing arrays containing coverages

        for r in samfile.fetch( until_eof = True ) : # for all records in the given BAM file
            if r.mapq > thres_mapq : # for record with mapq above the given threshold
                for start, end in r.blocks :
                    dict_refname_to_cov[ r.reference_name ][ start : end ] += 1 # increase coverage of aligned regions

    return dict_refname_to_cov
        