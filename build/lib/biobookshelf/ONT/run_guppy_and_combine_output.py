#!/usr/bin/env python
# 20210410

from biobookshelf.main import *
from biobookshelf import *
import argparse
import os, sys, getopt
from io import StringIO
import time
import math


def main( ) :
    parser = argparse.ArgumentParser( description = "Run Guppy basecaller on the nanopore sequencing datafiles in the given folder . This program has been developed by Hyunsu An (2021/04/09)." )
    parser.add_argument( "-d", "--dir_folder_nanopore_sequencing_data", help = "(Required) Run Guppy basecaller on the nanopore sequencing datafiles in the given folder" )
    parser.add_argument( "-b", "--flag_barcoding_was_used", help = "Set a flag indicating a barcoding kit was used during sequencing", action = 'store_true' )
    parser.add_argument( "-o", "--dir_folder_output_fastq", help = "(optional) copy fastq files in the guppy output folder to the given directory" )

    args = parser.parse_args( )
    if args.dir_folder_nanopore_sequencing_data is None :
        print( "required arguments are not given, exiting" )
        sys.exit( )

    # [input] parse arguments
    dir_folder_nanopore_sequencing_data = Program__Get_Absolute_Path_of_a_File( args.dir_folder_nanopore_sequencing_data )
    if dir_folder_nanopore_sequencing_data[ -1 ] != '/' : # add '/' to the end of the directory
        dir_folder_nanopore_sequencing_data += '/'
    flag_barcoding_was_used = args.flag_barcoding_was_used
    
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
        
    if args.dir_folder_output_fastq is not None : # if output folder of fastq files was given
        dir_folder_output_fastq = Program__Get_Absolute_Path_of_a_File( args.dir_folder_output_fastq )
        if dir_folder_output_fastq[ -1 ] != '/' : # add '/' to the end of the directory
            dir_folder_output_fastq += '/'
        for dir_file_fastq_gz in l_dir_file_fastq_gz : # for each output fastq file, copy the file to the given fastq folder
            shutil.copyfile( dir_file_fastq_gz, f"{dir_folder_output_fastq}{dir_file_fastq_gz.rsplit( '/', 1 )[ 1 ]}" )
        
if __name__ == "__main__" :
    main( )
    