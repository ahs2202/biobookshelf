from biobookshelf.main import *

def FASTA_Parse_UniProt_Header( line ) :
    """ # 2021-04-06 11:51:13 
    Parse uniprot header and return a dictionary containing field name and field values
    """
    if line[ 0 ] == '>' : # remove '>' character in front of the line if presents
        line = line[ 1 : ]
    l_split = list( r.replace( '_|_', ' = ' ) for r in line.strip( ).replace( ' = ', '_|_' ).split( '=' ) ) # some OS name contain ' = ' characters. mask ' = ' before spliting with '='. 
    l = [ ]
    for r in l_split[ : -1 ] :
        l.extend( r.rsplit( ' ', 1 ) )
    arr_r = np.array( [ 'UniProt_FASTA_Name' ] + l + l_split[ -1 : ], dtype = object )
    dict_field = dict( ( name_field, val ) for name_field, val in arr_r.reshape( ( int( len( arr_r ) / 2 ), 2 ) ) )
    dict_field[ 'header' ] = line # add original line
    return dict_field