#!/usr/bin/env python
# 20210401

from biobookshelf.main import *
from biobookshelf import *
import argparse
import os, sys, getopt
from io import StringIO
import time
import math


def main( ) :
    parser = argparse.ArgumentParser( description = "Analyze plasmid sequence with nanopore sequencing. This program has been developed by Hyunsu An (2021/04/01)." )
    parser.add_argument( "-r", "--dir_file_fasta_ref", help = "(Required) directory of a fasta file containing the reference sequence. (e.g. a fasta sequence from AddGene)" )
    parser.add_argument( "-i", "--dir_file_fastq", help = "(Required) directory of a fastq file from a nanopore sequencing" )
    parser.add_argument( "-o", "--dir_folder_output", help = "(Default: subdirectory of the folder containing the given fastq file) directory of output folder", default = 'default' )
    parser.add_argument( "-t", "--threads", help = "(Default: 10) Number of threads to use in the current compute node.", default = '10' )

    args = parser.parse_args( )
    if args.dir_file_fasta_ref is None or args.dir_file_fastq is None  :
        print( "required arguments are not given, exiting" )
        sys.exit( )

    # [input] parse arguments
    dir_file_fasta_ref = Program__Get_Absolute_Path_of_a_File( args.dir_file_fasta_ref )
    dir_file_fastq = Program__Get_Absolute_Path_of_a_File( args.dir_file_fastq )
    dir_folder_output = args.dir_folder_output
    n_threads = int( args.threads )

    # [input] output folder 
    # set default output folder 
    if dir_folder_output == 'default' : 
        dir_folder, name_file = dir_file_fastq.rsplit( '/', 1 )
        dir_folder += '/'
        if name_file.rsplit( '.', 1 )[ 1 ].lower( ) == 'gz' :
            name_file = name_file.rsplit( '.', 1 )[ 0 ]
        name_file = name_file.rsplit( '.', 1 )[ 0 ]
        dir_folder_output = f"{dir_folder}{name_file}/"

    dir_folder_output = Program__Get_Absolute_Path_of_a_File( dir_folder_output )
    if os.path.exists( dir_folder_output ) : # if output folder already exists, exit
        print( 'output folder already exists, exiting' )
        sys.exit( )
    else : # create an output folder
        os.makedirs( dir_folder_output )
    if dir_folder_output[ -1 ] != '/' : # add '/' at the end of the output folder
        dir_folder_output += '/'

    dir_file_index_ref = f"{dir_folder_output}index.ont.mmi"
    ONT.Minimap2_Index( dir_file_fasta_ref, dir_file_index_ref )
    dir_folder_minimap2_output = f"{dir_folder_output}minimap2/"
    ONT.Minimap2_Align( dir_file_fastq, dir_file_minimap2_index = dir_file_index_ref, dir_folder_minimap2_output = dir_folder_minimap2_output, n_threads = n_threads )

    dict_fasta = FASTA_Read( dir_file_fasta_ref )
    str_fasta_ref = dict_fasta[ list( dict_fasta )[ 0 ] ].upper( ) # read genome sequence of the reference (first sequence of the reference fasta file)
    shutil.copyfile( dir_file_fasta_ref, f"{dir_folder_minimap2_output}{dir_file_fasta_ref.rsplit( '/', 1 )[ 1 ]}" ) # copy reference sequence file to minimap2 output file (convenient minimap visualization

    # count number of each base for each position of the reference genome
    str_list_of_bases = "ATGC-"
    dict_base_to_index = dict( ( ( str_base, i ) for i, str_base in enumerate( str_list_of_bases ) ) )
    arr_count_base = np.zeros( ( 5, len( str_fasta_ref ) ) )

    l = [ ]
    with pysam.AlignmentFile( glob.glob( f"{dir_folder_minimap2_output}*.bam" )[ 0 ], 'rb' ) as samfile :
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
    df_summary.to_excel( f"{dir_folder_output}summary__base_frequency.xlsx" )

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
    df_substitution.to_excel( f"{dir_folder_output}summary__substitution.xlsx" )
    
if __name__ == "__main__" :
    main( )
    