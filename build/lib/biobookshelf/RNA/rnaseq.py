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
        
def STAR_Align( dir_folder_index, dir_file_read_1, dir_file_read_2 = None, dir_prefix_output = None, n_threads = 10, index_bam = True, star_option = None ) :
    """
    # 2021-07-05 21:34:33 
    Align reads with STAR aligner
    
    'dir_folder_index' : dir_folder containing STAR index
    'dir_file_read_1' : file containing read_1 (gzipped file supported)
    'dir_file_read_2' : file containing read_2 (gzipped file supported). read_1 file and read_2 file should be all gzipped or unzipped.
    'dir_prefix_output' : prefix of the STAR output files
    'index_bam' : index bam file using samtools
    'star_option' : additional command line option for star (arguments separated by spaces)
    
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
    if star_option is not None :
        l_args += star_option.strip( ).split( )
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
    
def Variant_Calling( r, dict_genome, set_mut_filtered = None, return_corrected_read_sequence = False, flag_ignore_indel = False, flag_return_as_string = True, flag_return_matched = False ) :
    ''' # 2021-08-27 20:46:12 
    perform variant calling of a single aligned read (pysam read object) using a given genome
    return a list of mutations and corrected read sequence
    'set_mut_filtered' : only consider mutations in the given 'set_mut_filtered'. only valid when 'return_corrected_read_sequence' = True
    'return_corrected_read_sequence' : return aligned read sequence after hard-clipping and filtering variants using 'set_mut_filtered'
    'flag_ignore_indel' : ignore insertion and deletion mutation calls
    'flag_return_as_string' : return a called mutation as a string in the following format: f"{reference_name}:{start}_{str_operation}_{mutation}". If set to False, return the called mutation as a tuple (reference_name, start, str_operation, mutation, quality)
    'flag_return_matched' : return 'id_mut' even when there is no mutation. the 'str_operation' for matched base is 'M'
    '''
    pos_read, pos_ref = 0, 0 # current position in the extracted reference sequence and the current read
    str_seq_ref = dict_genome[ r.reference_name ][ r.reference_start : r.reference_start + r.alen ] # retrieve a part of the reference sequence where the current read was aligned
    l_mut = [ ]
    
    if return_corrected_read_sequence : # correction mode
        str_seq_corrected_read = '' # corrected read sequence after hard clipping and filtering variants
        for n_bases, str_oper in NGS_CIGAR_Iterate_a_CIGAR_String( r.cigarstring ) :
            if str_oper == 'M' :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = r.seq[ pos_read ], str_seq_ref[ pos_ref ], r.query_qualities[ pos_read ]
                    str_base_read_corrected = str_base_ref # default corrected read base = ref
                    if str_base_read != str_base_ref :
                        id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_X_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'X', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        if set_mut_filtered is None or id_mut in set_mut_filtered :
                            l_mut.append( id_mut )
                            str_base_read_corrected = str_base_read
                    elif flag_return_matched :
                        id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_M_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                        l_mut.append( id_mut )
                    str_seq_corrected_read += str_base_read_corrected # add corrected read 
                    pos_ref += 1
                    pos_read += 1
            elif str_oper == 'N' :
                pos_ref += n_bases
            elif str_oper == 'S' :
                pos_read += n_bases
            elif str_oper == 'I' :
                id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_I_{r.seq[ pos_read : pos_read + n_bases ]}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'I', r.seq[ pos_read : pos_read + n_bases ], r.query_qualities[ pos_read : pos_read + n_bases ] ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and ( set_mut_filtered is None or id_mut in set_mut_filtered ) :
                    l_mut.append( id_mut ) 
                    str_seq_corrected_read += r.seq[ pos_read : pos_read + n_bases ] # add inserted sequence (if insertion is valid)
                pos_read += n_bases
            elif str_oper == 'D' :
                id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_D_{n_bases}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'D', n_bases, None ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                if not flag_ignore_indel and ( set_mut_filtered is None or id_mut in set_mut_filtered ) :
                    l_mut.append( id_mut ) 
                else :
                    str_seq_corrected_read += str_seq_ref[ pos_ref : pos_ref + n_bases ] # add deleted reference sequence (if deletion is invalid)
                pos_ref += n_bases
            elif str_oper == 'H' :
                pass
            elif str_oper == 'P' :
                pass
            elif str_oper == '=' :
                if flag_return_matched :
                    str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
                    id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_M_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append( id_mut )
                str_seq_corrected_read += r.seq[ pos_read : pos_read + n_bases ] # add read sequences
                pos_ref += n_bases
                pos_read += n_bases
            elif str_oper == 'X' :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = r.seq[ pos_read ], str_seq_ref[ pos_ref ], r.query_qualities[ pos_read ]
                    str_base_read_corrected = str_base_ref # default corrected read base = ref
                    id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_X_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'X', str_base_read, str_qual_read ) # compose id_mut, 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    if set_mut_filtered is None or id_mut in set_mut_filtered :
                        l_mut.append( id_mut ) 
                        str_base_read_corrected = str_base_read
                    str_seq_corrected_read += str_base_read_corrected # add corrected read 
                    pos_ref += 1
                    pos_read += 1
        return l_mut, str_seq_corrected_read
    
    else : # normal mutation calling mode
        for n_bases, str_oper in NGS_CIGAR_Iterate_a_CIGAR_String( r.cigarstring ) :
            if str_oper == 'M' :
                for _ in range( n_bases ) :
                    str_base_read, str_base_ref, str_qual_read = r.seq[ pos_read ], str_seq_ref[ pos_ref ], r.query_qualities[ pos_read ]
                    if str_base_read != str_base_ref :
                        l_mut.append( f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_X_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'X', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    elif flag_return_matched :
                        l_mut.append( f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_M_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'M', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
            elif str_oper == 'N' :
                pos_ref += n_bases
            elif str_oper == 'S' :
                pos_read += n_bases
            elif str_oper == 'I' :
                if not flag_ignore_indel :
                    l_mut.append( f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_I_{r.seq[ pos_read : pos_read + n_bases ]}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'I', r.seq[ pos_read : pos_read + n_bases ], r.query_qualities[ pos_read : pos_read + n_bases ] ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_read += n_bases
            elif str_oper == 'D' :
                if not flag_ignore_indel :
                    l_mut.append( f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_D_{n_bases}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'D', n_bases, None ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                pos_ref += n_bases
            elif str_oper == 'H' :
                pass
            elif str_oper == 'P' :
                pass
            elif str_oper == '=' :
                if flag_return_matched :
                    str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
                    id_mut = f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_M_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'M', str_base_read, str_qual_read ) # compose id_mut, convert 0-based coordinate to 1-based coordinate when calling a mutation as a string
                    l_mut.append( id_mut )
                pos_ref += n_bases
                pos_read += n_bases
            elif str_oper == 'X' :
                for _ in range( n_bases ) :
                    str_base_read, str_qual_read = r.seq[ pos_read ], r.query_qualities[ pos_read ]
                    l_mut.append( f"{r.reference_name}:{r.reference_start + 1 + pos_ref}_X_{str_base_read}" if flag_return_as_string else ( r.reference_name, r.reference_start + pos_ref, 'X', str_base_read, str_qual_read ) ) # 0-based coordinate to 1-based coordinate when returning id_mut as a string
                    pos_ref += 1
                    pos_read += 1
        return l_mut




