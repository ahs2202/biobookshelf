from biobookshelf.main import *

def Download_Data( dir_file, dir_remote, name_package ) :
    """
    Download large datafile of a PyPI package from a remote directory. If the file already exists, skip downloading and exit the funciton.
    download given file (relative directory within a package) from the remote directory (usually a Github directory, with sufficient number of data packs for Git-LFS pulls) to a current package
    """
    dir_file_download_flag = pkg_resources.resource_filename( name_package, f'{dir_file}.download_completed.flag' )
    if not os.path.exists( dir_file_download_flag ) :
        OS_Download( f'{dir_remote}{dir_file}', pkg_resources.resource_filename( name_package, dir_file ) )
        with open( dir_file_download_flag, 'w' ) as file :
            file.write( 'download completed at ' + TIME_GET_timestamp( ) )
        print( f"data file {dir_file} of the package {name_package} has been downloaded." )
        
def Gunzip_Data( dir_file, name_package ) :
    """
    Unzip large datafiles of a PyPI package often downloaded from a remote directory. If the gunzipped file already exists, skip unzipping and exit the funciton.
    """
    if dir_file.count( '.' ) == 0 or dir_file.rsplit( '.', 1 )[ 1 ].lower( ) != 'gz' : # check whether given file is a gzipped file
        print( f"[Gunzip_Data] given file {dir_file} seems to be not a gzipped file, exiting" )
        return -1
    dir_file_download_flag = pkg_resources.resource_filename( name_package, f'{dir_file}.gunzip_completed.flag' )
    if not os.path.exists( dir_file_download_flag ) : # if gunzip has not been done or not completed 
        
        dir_file_zipped = pkg_resources.resource_filename( name_package, dir_file ) # absolute directory of gzipped file
        dir_file_unzipped = dir_file_zipped.rsplit( '.', 1 )[ 0 ] # retrieve directory of an unzipped file
        if os.path.exists( dir_file_unzipped ) : # if incomplete unzipped file exists, remove the file
            os.remove( dir_file_unzipped )
        OS_Run( [ 'gunzip', dir_file_zipped ] )
        with open( dir_file_download_flag, 'w' ) as file :
            file.write( 'gunzip completed at ' + TIME_GET_timestamp( ) )
        print( f"data file {dir_file} of the package {name_package} has been gunzipped." )
        