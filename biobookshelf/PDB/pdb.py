import numpy as np
import pandas as pd
from io import StringIO 
import collections # for external function
from copy import deepcopy

# import internal functions
def PD_Select( df, deselect = False, ** dict_select ) :
    ''' Select and filter rows of df according to the given dict_select. If 'deselect' is set to True, deselect rows according to the given dict_select  Usage example : PANDAS_Select( df_meta_imid_ubi, dict(  Data_Type = [ 'Proteome', 'Ubi_Profiling' ], Value_Type = 'log2fc' ) ) '''
    for col, query in dict_select.items( ) :
        if type( df ) is pd.Series :
            data_values = df.index.values if col == 'index' else df.values # select values or indices of a given pd.Series
        elif type( df ) is pd.DataFrame : 
            if col not in df.columns.values and col != 'index' :
                print( "'{}' does not exist in columns of a given DataFrame".format( col ) )
                continue
            data_values = df.index.values if col == 'index' else df[ col ].values
        else :
            print( '[INVALID INPUT]: Inputs should be DataFrame or Series' )
            return -1
        if isinstance( query, ( list, tuple, np.ndarray, set ) ) : # if data to be selected is iterable
            query = set( query ) if isinstance( query, ( list, tuple, np.ndarray ) ) else query  # convert query into set
            df = df[ list( False if data_value in query else True for data_value in data_values ) ] if deselect else df[ list( True if data_value in query else False for data_value in data_values ) ]
        else :
            df = df[ data_values != query ] if deselect else df[ data_values == query ]
    return df

def LIST_COUNT( iterable, return_series = True, duplicate_filter = 2, dropna = True, sort_series_by_values = True, convert_tuple_to_string = False ) :
    ''' 
    # 20210224
    return a dictionary where key = each element in a given list, value = counts of the element in the list. if 'duplicate_filter' is not None, return entries that are duplicated 'duplicate_filter' times or more. '''
    if dropna and isinstance( iterable, pd.Series ) : iterable = iterable.dropna( ) # if dropna is set to 'True', dropn NaN values before counting
    if isinstance( next( iterable.__iter__( ) ), ( np.ndarray, list ) ) : iterable = list( map( tuple, iterable ) ) # if value is non-hashable list of numpy array, convert a value to a hashable format, tuple
    dict_counted = COUNTER( iterable )
    if convert_tuple_to_string : # if 'convert_tuple_to_string' is True and values in a given list are tuples, convert tuples into string
        dict_counted__tuple_converted_to_string = dict( )
        for key in dict_counted :
            value = dict_counted[ key ]
            if isinstance( key, ( tuple ) ) : dict_counted__tuple_converted_to_string[ ( '{}, ' * len( key ) )[ : -2 ].format( * key ) ] = value # convert tuple into string concatanated with ', '
            else : dict_counted__tuple_converted_to_string[ key ] = value
        dict_counted = dict_counted__tuple_converted_to_string
    if return_series :
        s_counted = pd.Series( dict_counted )
        if duplicate_filter is not None : s_counted = s_counted[ s_counted >= duplicate_filter ]
        if sort_series_by_values : s_counted = s_counted.sort_values( ascending = False )
        return s_counted
    else : return dict_counted

# In[7]:


df_meta_pdb_format = pd.read_csv( StringIO( '''Columns	Data	Justification	Data Type
1-6	ATOM_or_HETATM	.	character
7-11	Atom_serial_number	right	integer
13-16	Atom_name	left*	character
17	Alternate_location_indicator		character
18-20	Residue_name	right	character
22	Chain_identifier		character
23-26	Residue_sequence_number	right	integer
27	Code_for_insertions_of_residues		character
31-38	X_orthogonal_Å_coordinate	right	real (8.3)
39-46	Y_orthogonal_Å_coordinate	right	real (8.3)
47-54	Z_orthogonal_Å_coordinate	right	real (8.3)
55-60	Occupancy	right	real (6.2)
61-66	Temperature_factor	right	real (6.2)
73-76	Segment_identifier	left	character
77-78	Element_symbol	right	character
79-80	Charge	.	character''' ), sep = '\t' )


# In[3]:


df_amino_acids = pd.read_csv( StringIO( '''Amino_acid	Letter_code__3	Letter_code__1	Class	Polarity	Charge__at_pH_7_4	Hydropathy_index	Molar_absorptivity__Wavelength__λmax__nm_	Molar_absorptivity__Coefficient__ε__mM_1·cm_1_	Molecular_mass	Abundance_in_proteins___Percent___129_	Standard_genetic_coding__IUPAC_notation	MaxASA__Tien_et_al__2013__theor__	MaxASA__Tien_et_al__2013__emp__	MaxASA__Miller_et_al__1987	MaxASA__Rose_et_al__1985
Alanine	Ala	A	Aliphatic	Nonpolar	Neutral	1.8			89.094	8.76	GCN	129.0	121.0	113.0	118.1
Arginine	Arg	R	Basic	Basic polar	Positive	−4.5			174.203	5.78	MGR, CGY (coding codons can also be expressed by: CGN, AGR)	274.0	265.0	241.0	256.0
Asparagine	Asn	N	Amide	Polar	Neutral	−3.5			132.119	3.93	AAY	195.0	187.0	158.0	165.5
Aspartate	Asp	D	Acid	Acidic polar	Negative	−3.5			133.10399999999998	5.49	GAY	193.0	187.0	151.0	158.7
Cysteine	Cys	C	Sulfuric	Nonpolar	Neutral	2.5	250	0.3	121.154	1.38	UGY	167.0	148.0	140.0	146.1
Glutamate	Glu	E	Acid	Acidic polar	Negative	−3.5			147.131	6.32	GAR	223.0	214.0	183.0	186.2
Glutamine	Gln	Q	Amide	Polar	Neutral	−3.5			146.14600000000002	3.9	CAR	225.0	214.0	189.0	193.2
Glycine	Gly	G	Aliphatic	Nonpolar	Neutral	−0.4			75.067	7.03	GGN	104.0	97.0	85.0	88.1
Histidine	His	H	Basic aromatic	Basic polar	Positive, 10%, Neutral, 90%	−3.2	211	5.9	155.156	2.26	CAY	224.0	216.0	194.0	202.5
Isoleucine	Ile	I	Aliphatic	Nonpolar	Neutral	4.5			131.175	5.49	AUH	197.0	195.0	182.0	181.0
Leucine	Leu	L	Aliphatic	Nonpolar	Neutral	3.8			131.175	9.68	YUR, CUY (coding codons can also be expressed by: CUN, UUR)	201.0	191.0	180.0	193.1
Lysine	Lys	K	Basic	Basic polar	Positive	−3.9			146.189	5.19	AAR	236.0	230.0	211.0	225.8
Methionine	Met	M	Sulfuric	Nonpolar	Neutral	1.9			149.208	2.32	AUG	224.0	203.0	204.0	203.4
Phenylalanine	Phe	F	Aromatic	Nonpolar	Neutral	2.8	257, 206, 188	0.2, 9.3, 60.0	165.192	3.87	UUY	240.0	228.0	218.0	222.8
Proline	Pro	P	Cyclic	Nonpolar	Neutral	−1.6			115.132	5.02	CCN	159.0	154.0	143.0	146.8
Serine	Ser	S	Hydroxylic	Polar	Neutral	−0.8			105.09299999999999	7.14	UCN, AGY	155.0	143.0	122.0	129.8
Threonine	Thr	T	Hydroxylic	Polar	Neutral	−0.7			119.119	5.53	ACN	172.0	163.0	146.0	152.5
Tryptophan	Trp	W	Aromatic	Nonpolar	Neutral	−0.9	280, 219	5.6, 47.0	204.22799999999998	1.25	UGG	285.0	264.0	259.0	266.3
Tyrosine	Tyr	Y	Aromatic	Polar	Neutral	−1.3	274, 222, 193	1.4, 8.0, 48.0	181.19099999999997	2.91	UAY	263.0	255.0	229.0	236.8
Valine	Val	V	Aliphatic	Nonpolar	Neutral	4.2			117.148	6.73	GUN	174.0	165.0	160.0	164.5''' ), sep = '\t' )


# In[4]:


df_amino_acids


# In[1]:


"""def Parse_ATOM_or_HETATM_line( line ) :
    '''  parse a given line (should be ATOM or HETATM records) and return parsed data values as a numpy array (dtype = object)  '''
    try :
        arr_values = np.array( [ line[ 0 : 6 ].strip( ), int( line[ 6 : 11 ].strip( ) ), line[ 12 : 16 ].strip( ), line[ 16 ].strip( ), line[ 17 : 20 ].strip( ), line[ 21 ].strip( ), int( line[ 22 : 26 ].strip( ) ), line[ 26 ].strip( ), float( line[ 30 : 38 ].strip( ) ), float( line[ 39 : 46 ].strip( ) ), float( line[ 46 : 54 ].strip( ) ), float( line[ 54 : 60 ].strip( ) ), float( line[ 60 : 66 ].strip( ) ), line[ 72 : 76 ].strip( ), line[ 76 : 78 ].strip( ), line[ 78 : 80 ].strip( ) ], dtype = object )
    except : # for ZDOCK result (contain four wrongly placed columns containing real values)
        arr_values = np.array( [ line[ 0 : 6 ].strip( ), int( line[ 6 : 11 ].strip( ) ), line[ 12 : 16 ].strip( ), line[ 16 ].strip( ), line[ 17 : 20 ].strip( ), line[ 21 ].strip( ), int( line[ 22 : 26 ].strip( ) ), line[ 26 ].strip( ), float( line[ 30 : 38 ].strip( ) ), float( line[ 39 : 46 ].strip( ) ), float( line[ 46 : 54 ].strip( ) ), line[ 54 : 60 ].strip( ), line[ 60 : 66 ].strip( ), line[ 72 : 76 ].strip( ), line[ 76 : 78 ].strip( ), line[ 78 : 80 ].strip( ) ], dtype = object )
    return arr_values""" # 20191104
"""def Compose_ATOM_or_HETATM_line( arr_values ) :
    '''  compose a ATOM or HETATM line from the output of 'Parse_ATOM_or_HETATM_line' (it is essentially a reverse operation of 'Parse_ATOM_or_HETATM_line') '''
    return "{:<6s}{:>5d} {:<4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<3s}{:>3s}{:2s}".format( * arr_values )""" # 20191104


# In[4]:


def Parse_ATOM_or_HETATM_line( line ) :
    '''  # 2021-04-27 11:12:31 
    parse a given line (should be ATOM or HETATM records) and return parsed data values as a numpy array (dtype = object) Parse only until xyz coordinates (because later columns are somewhat variable across different docking servers)  '''
    # parse optional fields based on the length of the entry
    value = '' if len( line ) < 60 else line[ 54 : 60 ].strip( )
    float_occupancy = np.nan if len( value ) == 0 else float( value )
    value = '' if len( line ) < 66 else line[ 60 : 66 ].strip( )
    float_temperature_factor = np.nan if len( value ) == 0 else float( value )
    str_id_segment = '' if len( line ) < 76 else line[ 72 : 76 ].strip( )
    str_symbol_element = '' if len( line ) < 78 else line[ 76 : 78 ].strip( )
    str_change = '' if len( line ) < 80 else line[ 78 : 80 ].strip( )
    # parse required fields
    arr = np.array( [ line[ 0 : 6 ].strip( ), None if "*" in line[ 6 : 11 ] else int( line[ 6 : 11 ].strip( ) ), line[ 12 : 16 ].strip( ), line[ 16 ].strip( ), line[ 17 : 20 ].strip( ), line[ 21 ].strip( ), int( line[ 22 : 26 ].strip( ) ), line[ 26 ].strip( ), float( line[ 30 : 38 ].strip( ) ), float( line[ 39 : 46 ].strip( ) ), float( line[ 46 : 54 ].strip( ) ), float_occupancy, float_temperature_factor, str_id_segment, str_symbol_element, str_change ], dtype = object )
    # atom_number might contain non-integer characters ('*****') when the atom_number is > 99999 (e.g. pdb files from 'SWISS-MODEL repositories')
    return arr


# In[11]:


def Compose_ATOM_or_HETATM_line( arr_values ) :
    '''  # 2020-12-14 00:08:19 
    compose a ATOM or HETATM line from the output of 'Parse_ATOM_or_HETATM_line' (it is essentially a reverse operation of 'Parse_ATOM_or_HETATM_line') '''
    arr_values = deepcopy( arr_values )
    if len( df_meta_pdb_format ) > len( df_meta_pdb_format ) :
        arr_values = arr_values[ : len( df_meta_pdb_format ) ]
    arr_values[ 11 ] = "      " if isinstance( arr_values[ 11 ], ( str ) ) or np.isnan( arr_values[ 11 ] ) else "{:>6.2f}".format( arr_values[ 11 ] )
    arr_values[ 12 ] = "      " if isinstance( arr_values[ 12 ], ( str ) ) or np.isnan( arr_values[ 12 ] ) else "{:>6.2f}".format( arr_values[ 12 ] )
    
    atom_name = arr_values[ 2 ]
    if atom_name[ 0 ] in "ONCHSP" and len( atom_name ) < 4 : # based on description 'Atom names start with element symbols right-justified in columns 13-14 as permitted by the length of the name'
        return "{:<6s}{:>5d}  {:<3s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6s}{:>6s}      {:<4s}{:>2s}{:<2s}".format( * arr_values )
    else :
        return "{:<6s}{:>5d} {:<4s}{:1s}{:>3s} {:1s}{:>4d}{:1s}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6s}{:>6s}      {:<4s}{:>2s}{:<2s}".format( * arr_values )


# In[ ]:


def Read_Single_Module( dir_file, enable_automatic_correction = False, verbose = False ) : # 2020-07-07 15:59:00 
    '''   # 2021-04-27 11:25:12 
    Read PDB file with a single module (one chain A and one chain B, etc), and return dataframes of ATOM and HETATM records. If missing chain identifiers are detected, reassign chain identifiers based on chain transitions inferred by residue_numbers  '''
    if '/' in dir_file : 
        with open( dir_file, 'r' ) as file :
            l_lines = file.read( ).split( '\n' )
    else : # if a given 'dir_file' is not a directory (does not contain '/') read 'dir_file' as a string
        l_lines = dir_file.split( '\n' )
    l_lines = list( line for line in l_lines if 'HEADER' not in line and 'REMARK' not in line ) # remove header or remarks from the data
    n_terminations = len( list( line for line in l_lines if 'TER' in line[ : 3 ] ) ) # count the number of terminations
    # parse all lines containing ATOM or HETATM records
    l_arr = [ ]
    for line in l_lines :
        if 'ATOM' in line[ : 6 ] or 'HETATM' in line[ : 6 ] : # parse all lines containing ATOM or HETATM records
            arr = Parse_ATOM_or_HETATM_line( line )
            if arr[ 1 ] is None : # if atom_number contain invalid number 
                arr[ 1 ] = l_arr[ -1 ][ 1 ] + 1 # set atom_number as previous_atom_number + 1
            l_arr.append( arr )
    if verbose : print( "{} lines of ATOM or HETATM and {} terminations are found".format( str( len( l_arr ) ), n_terminations ) )
    df = pd.DataFrame( np.vstack( l_arr ), columns = df_meta_pdb_format.Data.values ).sort_values( [ 'Chain_identifier', 'Residue_sequence_number', 'Atom_serial_number' ], ignore_index = True )  # build dataframes of ATOM and HETATM records (except the five last data columns)
    arr_index_chain_transition = np.where( np.diff( df.Residue_sequence_number.values ) < 0 )[ 0 ] # retrive array of row indices where chain transiton exist (where new chain starts) (values in arr_index_chain_transition indicates locations of new chain starts)
    n_chains_in_pdb, n_chains_inferred = len( df.Chain_identifier.unique( ) ), len( arr_index_chain_transition ) + 1
    if n_chains_in_pdb < n_chains_inferred and enable_automatic_correction : # detect error in chain identifier assignment and try to correct it
        if verbose : print( "[Read_Single_Module] [ERROR:Chain Identifier Missing] {} chains were found in PDB files, but the redisue number trend says there are {} chain. Chain Identifiers will be re-assigned".format( str( n_chains_in_pdb ), str( n_chains_inferred ) ) )
        arr_chain_identifier = np.full( len( df ), 'A' )
        l_index_chain_transition = [ 0 ] + list( arr_index_chain_transition + 1 ) + [ len( df ) ]
        for index, str_new_chain_identifier in zip( np.arange( n_chains_inferred ), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' ) :
            arr_chain_identifier[ l_index_chain_transition[ index ] : l_index_chain_transition[ index + 1 ] ] = str_new_chain_identifier
        df.Chain_identifier = arr_chain_identifier
    return df


# In[ ]:


def Read_Multiple_Models( dir_file ) : # 2020-11-27 16:27:07 
    """ read pdb file of 'dir_file' containing multiple models (CABS-Dock output 'clus' and 'replica' files) and return l_df """
    l_df, l_l = list( ), list( )
    with open( dir_file ) as file : 
        while True :
            line = file.readline( )
            if len( line ) == 0 : break
            if line[ : 5 ] == 'MODEL' :
                name_model = line.strip( ).split( )[ 1 ]
                l_l = list( )
            elif 'HEADER' not in line and 'REMARK' not in line and ( 'ATOM' in line or 'HETATM' in line ) : # if current line contains ATOM or HETATM records
                l_l.append( Parse_ATOM_or_HETATM_line( line ) )
            elif line[ : 6 ] == 'ENDMDL' :
                df = pd.DataFrame( l_l, columns = df_meta_pdb_format.Data.values ).sort_values( [ 'Chain_identifier', 'Residue_sequence_number', 'Atom_serial_number' ], ignore_index = True )
                df[ 'MODEL' ] = name_model # assign id_model
                l_l = list( )
                l_df.append( df )
    return l_df


# In[ ]:


def Split_Multiple_Models( dir_file, dir_folder_output = None ) : # 2020-11-30 03:07:27 
    """ read pdb file of 'dir_file' containing multiple models (CABS-Dock output 'clus' and 'replica' files) and split the file into individual models without parsing and composing data values.
    'dir_folder_output' : output directory of split files. by default, it is where the pdb file is located """
    dir_folder, name_file = dir_file.rsplit( '/', 1 )
    if dir_folder_output is None : dir_folder_output = dir_folder # set default folder
    dir_prefix_output = dir_folder_output + name_file.rsplit( '.', 1 )[ 0 ]
    flag_writing = False # flag for writing
    with open( dir_file ) as file : 
        while True :
            line = file.readline( )
            if len( line ) == 0 : break
            if line[ : 5 ] == 'MODEL' :
                name_model = line.strip( ).split( )[ 1 ]
                newfile = open( dir_prefix_output + '.model_{name_model}.pdb'.format( name_model = name_model ), 'w' )
                flag_writing = True
                newfile.write( line )
            elif 'HEADER' not in line and 'REMARK' not in line and ( 'ATOM' in line or 'HETATM' in line ) and flag_writing : # if current line contains ATOM or HETATM records and currently a file is opened for writing
                newfile.write( line )
            elif line[ : 6 ] == 'ENDMDL' :
                newfile.close( )
                flag_writing = False


# In[4]:


def Write_Single_Module( df, dir_file, verbose = False ) : # 2020-07-07 15:59:03 
    '''   Add a 'TER' record automatically at the end of each chain, and write values in a given PDB dataframe (output of 'Read_Single_Module' function) as a PDB text file in the given directory 'dir_file'.  '''
    df = df.sort_values( [ 'Chain_identifier', 'Residue_sequence_number', 'Atom_name' ] ) # sort df before adding 'TER' records after the end of each chain and save the PDB dataframe as a PDB text file 
    df.Atom_serial_number = np.arange( 1, len( df ) + 1 ) # renumber 'Atom_serial_number'
    l_line = list( map( Compose_ATOM_or_HETATM_line, df.values ) ) # convert values in dataframe to string according to PDB text format.
    l_index_chain_transition = [ 0 ] + list( np.where( np.diff( df.Residue_sequence_number.values ) < 0 )[ 0 ] + 1 ) + [ len( df ) ] # retrive array of row indices where where new chain starts
    n_ter = len( l_index_chain_transition ) - 1
    if verbose : print( "{} chains detected. 'TER' record will be added".format( n_ter ) )
    l_line_with_ter = list( ) # add 'TER' records
    for index in np.arange( n_ter ) :
        l_line_with_ter.extend( l_line[ l_index_chain_transition[ index ] : l_index_chain_transition[ index + 1 ] ] + [ 'TER' ] )
    if '.pdb' not in dir_file.lower( ) : dir_file += '.pdb' # add pdb file extension if it does not exist in the given directory to the file to be written
    with open( dir_file, 'w' ) as file : file.write( '\n'.join( l_line_with_ter ) + '\nEND\n' )


# In[ ]:


def Retrive_Signatures_of_atoms( df, l_col = [ 'Atom_name', 'Residue_name', 'Chain_identifier', 'Residue_sequence_number' ] ) :
    '''  Input : output dataframe of 'Read_Single_Module'. Output : a set of signatures of atoms (tuple of values contained in the given array of columns, default = [ 'Atom_name', 'Residue_name', 'Chain_identifier', 'Residue_sequence_number' ] )  '''
    arr_atoms = df[ l_col ].values
    s_duplicated_signatures = L.Count( list( tuple( arr_values ) for arr_values in arr_atoms ) )
    if len( s_duplicated_signatures ) > 0 :
        print( '{} duplicated signatures of atom were found'.format( str( len( s_duplicated_signatures ) ) ) )
    return set( tuple( arr_values ) for arr_values in arr_atoms )


# In[ ]:


def Identify_Protein_and_Assign_Chain_identifier( df, ** dict_identifiable_chains ) :
    '''  Input : output dataframe of 'Read_Single_Module' as 'df' and keyworded arguments for desired chain identifier for each identifiable protein (example input: A = { 'start' : { 'ILE' : 47 }, 323 : 'CYS', 326 : 'CYS', 391 : 'CYS', 394 : 'CYS' } )
    Output : output dataframe of 'Read_Single_Module', but with modified chain identifiers for each identified proteins  '''
    l_df_a_chain = list( )
    for char_chain_identifier in df.Chain_identifier.unique( ) : # resolve merged multiper proteins (for example ACAT1 dimer)
        df_a_chain = PD.Select( df, Chain_identifier = char_chain_identifier )
        arr_residue_num = df_a_chain.Residue_sequence_number.unique( )
        n_residues = len( arr_residue_num )
        if n_residues == 391 * 2 : # for HawkDock ACAT1 dimer (their residue numbers are renumbered and needed to be separated)
            print( '[Identify_Protein_and_Assign_Chain_identifier] ACAT1 Dimer Splited' )
            df_a_chain.loc[ df_a_chain.Residue_sequence_number > 391, 'Chain_identifier' ] = 'Splited'
        l_df_a_chain.append( df_a_chain )
    df = pd.concat( l_df_a_chain )
    l_df_a_chain = list( )
    for char_chain_identifier in df.Chain_identifier.unique( ) :    
        df_a_chain = PD.Select( df, Chain_identifier = char_chain_identifier )
        for char_desired_chain_identifier, dict_sequence_signature in dict_identifiable_chains.items( ) :
            dict_sequence_signature = deepcopy( dict_sequence_signature )
            str_start_residue, int_start_residue_num = list( dict_sequence_signature.pop( 'start' ).items( ) )[ 0 ] # retrive given start residue number and name for the unique chain profile ('dict_sequence_signature') 
            int_min_residue_num = df_a_chain.Residue_sequence_number.min( )
            if int_min_residue_num != int_start_residue_num and PD.Select( df_a_chain, Residue_sequence_number = int_min_residue_num ).Residue_name.values[ 0 ] == str_start_residue : # if the starting residue names are the same but the starting residue numbers are different, re-number residue numbers of the chain
                print( "\t[Identify_Protein_and_Assign_Chain_identifier] [NOTE] the given start residue name and start residue name in the chain is the same while start residue numbers are different, and residue numbers will be re-numbered" )
                df_a_chain.Residue_sequence_number = df_a_chain.Residue_sequence_number + ( int_start_residue_num - int_min_residue_num )
            if np.sum( list( PD.Select( df_a_chain, Residue_sequence_number = int_residue_number ).Residue_name.values[ 0 ] != str_residue_name if len( PD.Select( df_a_chain, Residue_sequence_number = int_residue_number ) ) > 0 else True for int_residue_number, str_residue_name in dict_sequence_signature.items( ) ) ) == 0 :
                df_a_chain.Chain_identifier = char_desired_chain_identifier
                l_df_a_chain.append( df_a_chain )
                dict_identifiable_chains.pop( char_desired_chain_identifier )
                break
    return pd.concat( l_df_a_chain )


# In[ ]:


def Clean_Minimal( dir_file_pdb, dir_file_pdb_clean, l_chain_id = None, remove_hetatm = True ) : # 2020-11-30 01:57:22 
    ''' Clean a pdb (from RCSB PDB) in a minimal fashion
    'l_chain_id' : list of chain_ids to retain. by default, retain all chain_ids
    'remove_hetatm' : remove hetatm lines '''
    if l_chain_id is not None : set_chain_id = set( l_chain_id )
    with open( dir_file_pdb_clean, 'w' ) as newfile :
        with open( dir_file_pdb ) as file :
            while True :
                line = file.readline( )
                if len( line ) == 0 : break
                if remove_hetatm :
                    if line[ : 6 ] == 'HETATM' : continue
                if l_chain_id is not None :
                    if line[ 21 ] not in set_chain_id : continue
                newfile.write( line )


# In[ ]:


def Dissociate( df, * l_l_chains, float_distance_for_dissociation = 1000 ) :
    '''  # 2020-12-21 23:02:46 
    Move apart the parts defined by a list of lists of chain_ids 'l_l_chains' by 'float_distance_for_dissociation' from the center of the system.
    float_distance_for_dissociation = 1000 # distance used for dissociation of each 'part'
    '''
    l_col_coord = [ 'X_orthogonal_Å_coordinate', 'Y_orthogonal_Å_coordinate', 'Z_orthogonal_Å_coordinate' ]
    df = deepcopy( df )
    
    coord_center = df[ l_col_coord ].values.mean( axis = 0 ) # retrieve coordinates of center (average coordinates of all atoms)
    l_df = list( )
    for l_chains in l_l_chains :
        df_part = PD.Select( df, Chain_identifier = l_chains ) # subset pdb with the current list of chain_ids
        coord_center_part = df_part[ l_col_coord ].values.mean( axis = 0 ) # retrieve coordinates of center of the current part (average coordinates of all atoms)
        vector_center_part_from_center_system = coord_center_part - coord_center
        unit_vector_center_part_from_center_system = vector_center_part_from_center_system / sum( vector_center_part_from_center_system ** 2 ) ** 0.5 # retrieve unit vector of the centor of the part relative to the center of the system
        vector_move = unit_vector_center_part_from_center_system * float_distance_for_dissociation
        df_part.loc[ :, l_col_coord ] = df_part[ l_col_coord ].values + vector_move # move the part by the unit vector times 'float_distance_for_dissociation'
        l_df.append( df_part )
    df_moved = pd.concat( l_df ) # combine moved parts
    return df_moved


# # Specific Functions

# In[ ]:


def CRBN_Modify_Residue_Name( df ) :
    '''   Modify PDB file containing CRBN so that its Zn coordinating site can be appropriatly simulated. If PDB records were already modified, modify them back to standard residue names.   '''
    df = deepcopy( df )
    df = df[ df.Atom_name != 'ZN' ] # remove previous Zn from the structure before processing
    arr_mask_crbn = df.Chain_identifier.values == 'A'
    arr_residue_num = df.Residue_sequence_number.values
    if len( PD.Select( df, Residue_name = 'CYZ' ) ) == 0 : # if a given PDB dataframe's CRBN's residues were already modified.
        df.loc[ arr_mask_crbn & ( ( arr_residue_num == 323 ) | ( arr_residue_num == 326 ) | ( arr_residue_num == 391 ) ), 'Residue_name' ] = 'CY1'
        df.loc[ arr_mask_crbn & ( arr_residue_num == 394 ), 'Residue_name' ] = 'CYZ'
        df.loc[ -1 ] = np.array( [ 'ATOM', 3280, 'ZN', '', 'CYZ', 'A', 394, '' ] + list( PD.Select( df, Residue_name = [ 'CY1', 'CYZ' ], Atom_name = 'SG' )[ [ 'X_orthogonal_Å_coordinate', 'Y_orthogonal_Å_coordinate', 'Z_orthogonal_Å_coordinate' ] ].mean( axis = 0 ) ), dtype = object ) # add discarded Zn ion in the PDB dataframe (average of xyz coordinats of SG of the four cysteines)
    else :
        df.loc[ arr_mask_crbn & ( ( arr_residue_num == 323 ) | ( arr_residue_num == 326 ) | ( arr_residue_num == 391 ) | ( arr_residue_num == 394 ) ), 'Residue_name' ] = 'CYS'
    return df.sort_values( [ 'Chain_identifier', 'Residue_sequence_number', 'Atom_name' ] ).reset_index( drop = True )

