import biobookshelf.STR as STR
import biobookshelf.MAP as MAP

# basic python packages
import sys # to measure the size of an object
import os # Import library to interact with operating system
import subprocess # for running commands 
from subprocess import Popen, PIPE # for reading stdout
import glob # for getting file names
import inspect # module to access code content of during function call
import datetime # to retrive current time
import gzip # to handle gzip file
import base64 # for converting binary to text data (web application)
from io import StringIO # for converting a string to a file-like stream
import json # to read and write JSON file
import shutil # for copying file
import uuid # for universal identifier (same length with md5 hash) uuid.uuid4().hex
import gc # for explicit garbage collection
import time # for sleep function
from cycler import cycler # for cycling matplotlib color palette
from mpl_toolkits.mplot3d import Axes3D # module for 3D plotting
import numpy as np
import pandas as pd # read tab delimited file
import csv # for saving dataframe without extra "" quotations
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors # for normalization
from matplotlib import cm # to map scalar values to color using colormap
from matplotlib.collections import BrokenBarHCollection # for chromosome plotting
import math
import re # older version of regular expression module
# module for generating pairs for iteration
from itertools import combinations
import itertools
# import module for displaying Dataframe
# for adjusting jupyter notebook cell width
from IPython.core.display import display, HTML
# import module for copying object
from copy import deepcopy # deepcopy
from copy import copy # shallow_copy
import collections # count elements # usage : dict( collections.Counter( b))
from collections import defaultdict
import ast# usage : ast.literal_eval( string) # to convert string representation of list to list
import pickle # module for saving dictionary and dataframe (very fast when saving as binary!)
import traceback # printing traceback of exceptions
import mmap # for random access of file
from multiprocessing import Pool, get_context, set_start_method # for multiple processing  # with get_context("spawn").Pool() as pool:
import heapq # merge sorting files.
import contextlib # for use with heapq
import shlex, subprocess # for running multiple shell commands in python
## defining short cut for modules
np_str = np.core.defchararray
import importlib # for reloading custom modules after modifications
import requests # for retriving HTML documents
from ftplib import FTP # for interacting with ftp server
import urllib.request # to retrive html document from the internet
from xml.parsers.expat import ExpatError

# python modules that requires independent anaconda install
# import module for interactive plotting (with hover tool)
from bokeh.plotting import figure, output_file, show, output_notebook
# import module for opening URL by clicking entries
from bokeh.models import OpenURL, TapTool, ColumnDataSource, CustomJS, Rect, LabelSet, Label, Range1d, BoxAnnotation, LogColorMapper, LogTicker, ColorBar, Whisker
from bokeh.layouts import row, gridplot, widgetbox
from bokeh.transform import jitter
from bokeh.io import reset_output


# from bioservices.kegg import KEGG # module for KEGG REST service
# from Bio import SeqIO

import pysam # to read SAM and BAM file


import xmltodict # read xml as an ordered dictionary
# HTML
from bs4 import BeautifulSoup # for parsing HTML





# for plotly python
import plotly.express as px
import plotly 
import plotly.graph_objects as go

import plotnine as p9
import seaborn as sns


# single cell RNA-Seq data in Python
import scanpy as sc 
import scanpy.external as sce
import anndata

## binary arrays
from bitarray import bitarray
## for interval overlap searching
import intervaltree

# import statistical model for adjusting p-values
from statsmodels.stats.multitest import multipletests
import regex # regular expression modules for searching substrings in a string

# modules for calculating Pearson correlations 
import scipy
from scipy.stats import fisher_exact # usage : For GESA p_value, stats.fisher_exact([[80, 2000], [80, 200]])
from scipy.stats import mstats # import module for masked arrays
from scipy import stats
# from scipy import optimize as so # import a module for methods that can find minumum or maximum of functions 
import scipy # functions for hierarchial clusteringd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist # c, coph_dists = cophenet(Z, pdist(X)) #compares (correlates) the actual pairwise distances of all your samples to those implied by the hierarchical clustering.

import umap
from sklearn.decomposition import PCA, FactorAnalysis # import modules for multi-dimentional data visualization (t-SNE)
from sklearn.manifold import TSNE 
from sklearn.cluster import AgglomerativeClustering, KMeans, DBSCAN # # modules for clustering



def Wide( int_percent_html_code_cell_width = 95 ) :
    """ 
    # 20210224
    widen jupyter notebook cell
    """
    display(HTML("<style>.container { width:100% !important; }</style>"))
    display(HTML("<style>div.cell{width:" + str( int_percent_html_code_cell_width ) + "%;margin-left:" + str ( 100 - int_percent_html_code_cell_width ) + "%;margin-right:auto;}</style>" ) )

def Jupyter_Notebook_Extension__Build_Snippet( str_code, str_name, indentation = 2 ) :
    """ Print string that can be added to Snippet_menu Jupyter notebook extension """
    str_indent = '    ' * int( indentation )
    str_indent_code = '    ' * int( indentation + 2 )
    str_code_formatted = '["' + str_code.strip( ).replace( '"', r'\"' ).replace( "'", r'\"' ).replace( '\n', '",\n{str_indent_code}"'.format( str_indent_code = str_indent_code ) ) + '"]'
    str_json_block = str_indent + "{\n" +  '{str_indent}    "name" : "{str_name}",\n{str_indent}    "snippet" : {str_code_formatted}\n'.format( str_indent = str_indent, str_name = str_name, str_code_formatted = str_code_formatted ) + str_indent + "},\n"
    print( str_json_block )





# ### Define Reference Gene_Set for Normalization

# In[ ]:


known_reference_genes = [ 60, 203068, 11315, 59, 2597, 6161, 1915, 3251, 5430, 7846 ] # List_Gene__2__List_Gene_ID( ['ACTB', 'TUBB', 'PARK7', 'ACTA2', 'GAPDH', 'RPL32', 'EF1A', 'HPRT', 'POLR2A', 'TUBA1A' ] )
known_reference_genes_including_predicted = [ 60, 203068, 11315, 59, 2597, 6161, 1915, 3251, 5430, 7846, 25912, 27243, 56851, 2821, 5690, 5692, 7879, 7905, 6634,  7415, 51699 ]


# ### Useful Functions in applications 

# In[ ]:


def Program__Get_Absolute_Path_of_a_File( dir_file ) :
    ''' return an absolute path of a directory of a file '''
    if dir_file[ 0 ] == '/' : return dir_file # if dir_file is an absolute path, return dir_file without modification
    current_working_directory = os.getcwd( ) + '/'
    if dir_file[ : 2 ] == './' : return current_working_directory + dir_file[ 2 : ]
    else : return current_working_directory + dir_file


# In[ ]:


def UUID( ) :
    ''' return a 128bit universal unique identifier '''
    return uuid.uuid4( ).hex


# ### Functions common for general objects

# In[ ]:


def OBJECT_Get_size(obj, seen = None ):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([ OBJECT_Get_size(v, seen) for v in obj.values()])
        size += sum([ OBJECT_Get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += OBJECT_Get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([ OBJECT_Get_size(i, seen) for i in obj])
    return size


# ### Global functions

# In[ ]:


def GLOBALS_Find_object( query = None, case_specific = False ) : # 20190903
    ''' Search global object names with a given query. By default, return the list of all objects '''
    arr_object_names = np.array( list( globals( ).keys( ) ), dtype = object )
    if not case_specific : # perform case-insensitive search
        _, mask_matched = Search_list_of_strings( list( map( str.upper, arr_object_names ) ), query = query.upper( ), return_mask_matched = True ) if query is not None else l_object_names
        return list( arr_object_names[ mask_matched ] )
    return list( Search_list_of_strings( arr_object_names, query = query ) if query is not None else l_object_names )


# ### Functions for maintainance of Python workspace

# In[ ]:


def TIME_GET_timestamp( ) :
    '''  Get timestamp of current time in "%Y%m%d_%H%M" format  '''
    return datetime.datetime.now( ).strftime("%Y%m%d_%H%M")


# In[ ]:


############################ define objects that are required for some optimized methods #################################
arr_halfs = np.ones( 2 ) / 2
an_intercept_vector = np.array( [ [ 0.0 ], [ 0.0 ] ] ) # an empty 2D vector then will store the value of an intercept and become an intercept vector # for LINREGRESS
bbox_props_white = dict( boxstyle = "round", fc = "w", ec="0.5", alpha = 0.9 )  # define box properties that will be used in a function below for annotation

# ### General Python Functions

# In[ ]:


def _Get_object_name( data_object ) :
    ''' return the name of data_object '''
    return DATASET_retrive_name_of_data_object_of_outer_frame( '_Get_object_name' )


# In[ ]:


def CYCLER( l_entry ) :
    ''' cycle through given list (or a general iterable) '''
    n_entries = len( l_entry )
    index = 0
    while True :
        yield l_entry[ index ]
        index += 1
        if index == n_entries : index = 0


# ### Functions for handling integers

# In[ ]:


def INTEGER_Spread( n_values ) : # 2020-08-07 01:34:56 
    ''' spread 'n_values' number of integers ( 0 ~ n_values - 1 ) so that early integers are far apart from each other while later integers are adjacent to each other.
    Useful for assigning color on a 'jet' or 'rainbow' colormap '''
    step = n_values - 1
    index_cycle = 1
    int_encoding = 0 
    l_encoding = [ int_encoding ]
    while len( l_encoding ) < n_values : # until all integers are assigned
        int_encoding += step
        if int_encoding >= n_values : # detect the end of cycle, and start new cycle
            int_encoding -= n_values
            index_cycle += 1
            step = max( 1, int( n_values / math.pi ** ( index_cycle / 2 ) ) ) # use irrational number to define step for each cycle to spread integers more evenly # minimum step = 1INTEGER_Spread
            continue
        elif int_encoding in l_encoding : continue
        l_encoding.append( int_encoding )
    return l_encoding


# In[ ]:


def INT_Get_Ranges_from_List_of_Integers( l, flag_sorted = True ) :
    '''  # 2020-12-12 21:16:28 
    convert list of integers into list of ranges of integers
    '''
    if len( l ) == 0 : # return empty list of ranges if a given list is empty 
        return list( )
    if not flag_sorted : # sort the list if 'flag_sorted' is not True
        l = np.sort( l )
    diff_l = np.diff( l )
    l_index_break_points = [ -1 ] + list( np.where( diff_l != 1 )[ 0 ] ) + [ len( diff_l ) ] # identify start positions of break points
    l_range = list( )
    for i in range( len( l_index_break_points ) - 1 ) :
        l_range.append( [ l[ l_index_break_points[ i ] + 1 ], l[ l_index_break_points[ i + 1 ] ] ] )
    return l_range


# ### Functions for working with base systems.

# In[ ]:


def BASE( number, base, n_numbers = None ) : # 2020-07-25 16:08:23 
    ''' return a list of integers representing the number in the given base system.
    if 'n_numbers' is given and length of the output list is smaller than 'n_numbers', pad zeros and return a list of integers with the length of 'n_numbers' '''
    if isinstance( number, float ) : number = int( number )
    l_value = list( )
    quotient = number
    while True :
        quotient, remainder = divmod( quotient, base )
        l_value.append( remainder )
        if quotient == 0 : break
    if n_numbers is not None and len( l_value ) < n_numbers : # if 'n_numbers' is given and length of the output list is smaller than 'n_numbers', pad zeros and return a list of integers with the length of 'n_numbers'
        for index in range( n_numbers - len( l_value ) ) : l_value.append( 0 )
    return l_value[ : : -1 ]


# ### Functions for encoding data in a string/decoding data from a string

# In[ ]:


def Encode_List_of_Strings( arr, chr_separator = ';', chr_representing_repeated_string = '=' ) : # 2020-08-08 23:47:49 
    ''' Encode list of strings into a string. 
    'chr_representing_repeated_string': Reduce the size of string by representing repeating non-empty, valid string with single 'chr_representing_repeated_string' (default '=') character. To disable this behavior, set 'chr_representing_repeated_string' to None '''
    entry_previous = None
    l_value = list( )
    for entry in arr :
        if isinstance( entry, float ) : l_value.append( '' )
        elif chr_representing_repeated_string is not None and entry_previous == entry : l_value.append( chr_representing_repeated_string )
        else : 
            l_value.append( entry )
            entry_previous = entry
    return chr_separator.join( l_value )

def Decode_List_of_Strings( str_encoded, chr_separator = ';', chr_representing_repeated_string = '=' ) : # 2020-08-12 22:16:00 
    ''' Decode a string into a list of strings. Reduce sthe size of string by representing repeating non-empty, valid string with single 'chr_representing_repeated_string' (default '=') character '''
    entry_previous = None
    l_value = list( )
    for entry in str_encoded.split( chr_separator ) :
        if len( entry ) == 0 : l_value.append( np.nan )
        elif entry == chr_representing_repeated_string : l_value.append( entry_previous )
        else :
            l_value.append( entry )
            entry_previous = entry
    return np.array( l_value, dtype = object )


# ### Functions for workig with ASCII encodings

# In[ ]:


def ASCII_Decode( l_seq, ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 ) : # 2020-07-26 00:26:22 
    ''' Decode dictionary or list of strings containing encoded float values with one or more ascii character. (two ascii characters can encode float value with ~8000 levels)
    ascii_min: (integer) minimum = 33
    ascii_max: (integer) maximum = 126. should be larger than 'ascii_min' + 1. The last value will be used to interpret invalid values '''
    if ascii_max <= 125 : n_char = 1 # if not all ascii characters are used, use only one ascii character to encode/decode data
    elif n_char > 1 : ascii_max = 126
    base = ascii_max - ascii_min + 1 # base including ascii characters in l_ascii_to_exclude will be excluded 
    if l_ascii_to_exclude is not None : set_ascii_to_exclude = set( ascii_to_exclude for ascii_to_exclude in l_ascii_to_exclude if ascii_to_exclude - ascii_min >= 0 and ascii_to_exclude - ascii_min < base ) # retrieve set of ascii code to exclude in the encoding
    n_levels = ( base - len( set_ascii_to_exclude ) ) ** n_char - 1 if n_char > 1 else base - len( set_ascii_to_exclude ) - 1 # last level is for interpreting invalid value # ascii characters in l_ascii_to_exclude will be excluded in when interpreting ascii characters as levels
    step = ( value_max - value_min ) / n_levels # amount of an increment for each level
    arr_ascii_to_value = np.full( base ** n_char, np.nan ) # any ascii character outside the range of 'ascii_min' and 'ascii_max' and in 'set_ascii_to_exclude' will be considered as an invalid value (np.nan)
    value = value_min
    for index in range( len( arr_ascii_to_value ) - 1 ) : # last level is for interpreting invalid value (np.nan)
        set_ascii_used_for_representation = set( number + ascii_min for number in BASE( index, base, n_char ) )
        if sum( ascii_to_exclude in set_ascii_used_for_representation for ascii_to_exclude in set_ascii_to_exclude ) > 0 : continue
        arr_ascii_to_value[ index ] = value
        value += step
    bool_flag_input_is_dict = isinstance( l_seq, dict ) # booL_flag indicating whether input is given as a dictionary datatype
    l_arr_value = dict( ) if bool_flag_input_is_dict else list( )
    for seq in l_seq : # for each sequence
        if bool_flag_input_is_dict : key, seq = seq, l_seq[ seq ] # if dictionary has been given, retrieve key and sequence accordingly.
        len_seq = len( seq )
        arr_value = np.full( int( len_seq / n_char ), np.nan, dtype = float )
        for index in range( 0, len_seq, n_char ) :
            str_ascii = seq[ index : index + n_char ]
            if n_char > 1 :
                int_ascii_encoding = 0
                for char_ascii in str_ascii[ : -1 ] :
                    int_ascii_encoding += ord( char_ascii ) - ascii_min
                    int_ascii_encoding *= base
                int_ascii_encoding += ord( str_ascii[ -1 ] ) - ascii_min
            else : int_ascii_encoding = ord( str_ascii ) - ascii_min
            arr_value[ int( index / n_char ) ] = arr_ascii_to_value[ int_ascii_encoding ]
        if bool_flag_input_is_dict : l_arr_value[ key ] = arr_value
        else : l_arr_value.append( arr_value )
    return l_arr_value


# In[ ]:


def ASCII_Encode( l_l_value, ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 ) : # 2020-07-25 18:09:07 
    ''' Encode list of 'list of values' of float datatype with one or more ascii character, and return list of strings (ascii format) (each list of values in the given input will be encoded into a single string). (two ascii characters can encode float value with ~8000 levels)
    ascii_min: (integer) minimum = 33
    ascii_max: (integer) maximum = 126. should be larger than 'ascii_min' + 1. The last value will be used to interpret invalid values
    value_max: when encoded values are interger or floating number with regular intervalues, correct value_max should be the maximum integer value + step (1 for integer and amount of the regular interval for the float number), since values are floored down to an integer during the encoding process '''
    if ascii_max <= 125 : n_char = 1 # if not all ascii characters are used, use only one ascii character to encode/decode data
    elif n_char > 1 : ascii_max = 126
    base = ascii_max - ascii_min + 1 # base including ascii characters in l_ascii_to_exclude will be excluded 
    if l_ascii_to_exclude is not None : set_ascii_to_exclude = set( ascii_to_exclude for ascii_to_exclude in l_ascii_to_exclude if ascii_to_exclude - ascii_min >= 0 and ascii_to_exclude - ascii_min < base ) # retrieve set of ascii code to exclude in the encoding
    n_levels = ( base - len( set_ascii_to_exclude ) ) ** n_char - 1 if n_char > 1 else base - len( set_ascii_to_exclude ) - 1 # last level is for interpreting invalid value # ascii characters in l_ascii_to_exclude will be excluded in when interpreting ascii characters as levels
    str_ascii_encoding_representing_invalid_value = chr( ascii_max ) * n_char # an ascii_encoding string representing an invalid value
    arr_level_to_ascii_encoding = np.full( n_levels + 1, str_ascii_encoding_representing_invalid_value ) # initialize 'arr_level_to_ascii_encoding' filled with ascii_encoding representing an invalid value
    l_ascii = list( chr( value ) for value in range( ascii_min, ascii_max + 1, 1 ) if value not in set_ascii_to_exclude ) # list of valid ascii characters excluding 'set_ascii_to_exclude'
    for level in range( n_levels ) : # encode each level with available ascii characters (valid ascii characters, 'l_ascii')
        if n_char == 1 : arr_level_to_ascii_encoding[ level ] = l_ascii[ level ]
        else : arr_level_to_ascii_encoding[ level ] = ''.join( list( l_ascii[ number ] for number in BASE( level, len( l_ascii ), n_char ) ) ) # encode each level with available ascii characters (valid ascii characters, 'l_ascii')
    arr_level_to_ascii_encoding[ - 1 ] = arr_level_to_ascii_encoding[ - 2 ] # max + 1 level encode the same level as the max level
    l_seq = dict( ( l_value, ''.join( list( str_ascii_encoding_representing_invalid_value if np.isnan( value ) else arr_level_to_ascii_encoding[ int( ( value - value_min ) / ( value_max - value_min ) * n_levels ) ] for value in l_l_value[ l_value ] ) ) ) for l_value in l_l_value ) if isinstance( l_l_value, dict ) else list( ''.join( list( str_ascii_encoding_representing_invalid_value if np.isnan( value ) else arr_level_to_ascii_encoding[ int( ( value - value_min ) / ( value_max - value_min ) * n_levels ) ] for value in l_value ) ) for l_value in l_l_value ) # prepare the output according to the input datatype 
    return l_seq


# ### Functions for parsing file

# In[ ]:


def PARSE_Empty_Space_Delimited_with_Quotations( line ) : # 2020-07-17 01:04:25 
    """ parse empty-space (' ') delimited line. return only non-empty entries """
    char_quotation_mark = None # store quotation mark used in the quotation
    bool_flag_the_end_of_quotation = False # only true right after the end of a quotation
    str_value = ''
    l_value = list( )
    for character in line.strip( ) :
        if char_quotation_mark is None : # if quotation has not been encountered
            if character == ' ' : 
                l_value.append( str_value )
                str_value = ''
            elif character == "'" or character == '"' : 
                char_quotation_mark = character # if quotation mark is encountered
                str_value += character
            else : str_value += character
        else : # inside quotation
            if bool_flag_the_end_of_quotation : # if flag indicating the potential end of quotation has been raised
                if character == ' ' : # confirm the end of quotation 
                    bool_flag_the_end_of_quotation = False
                    l_value.append( str_value )
                    str_value = ''
                    char_quotation_mark = None
                elif character == char_quotation_mark : str_value += character # false positive (the flag still raised)
                else :
                    bool_flag_the_end_of_quotation = False # false positive (the flag is withdrawn)
                    str_value += character
            elif character == char_quotation_mark : # raise a flag indicating the potential end of quotation
                bool_flag_the_end_of_quotation = True
                str_value += character
            else : str_value += character
    l_value.append( str_value )
    return list( value for value in l_value if len( value ) > 0 ) # return only non-empty entries


# ### Functions for dictionary

# In[ ]:


def DICTIONARY_convert_a_dict_into_key_and_value_lists( a_dict ) :
    ''' return list_keys, list_values '''
    if type( a_dict ) is dict : # check a given object is a dictionary
        list_keys = list( )
        list_values = list( )
        for key, value in a_dict.items( ) :
            list_keys.append( key )
            list_values.append( value )
        return list_keys, list_values
    else : # if a_dict is not a dictionary, return an error value
        return -1 


# In[ ]:


def DICTIONARY_Build_from_arr( arr, order_index_entry = True ) :
    if order_index_entry :
        return dict( ( index, entry ) for entry, index in zip( arr, np.arange( len( arr ) ) ) ) 
    else :
        return dict( ( entry, index ) for entry, index in zip( arr, np.arange( len( arr ) ) ) ) 


# In[ ]:


def DICTIONARY_Find_Max( dict_value ) : # 2020-07-29 23:55:58 
    ''' find key value with the maximum value in a given dictionary, and return 'key_max', 'value_max' '''
    if len( dict_value ) == 0 : return None, None # if an empty dictionary is given, return None
    key_max = next( dict_value.__iter__( ) ) # retrieve the a key and value from the dictionary as an initial key-value pair with 'max' value
    value_max = dict_value[ key_max ]
    for key in dict_value :
        value = dict_value[ key ]
        if value_max < value :
            key_max = key
            value_max = value
    return key_max, value_max

def DICTIONARY_Find_Min( dict_value ) : 
    ''' 
    # 2020-12-06 18:33:24 
    find key value with the maximum value in a given dictionary, and return 'key_min', 'value_min' '''
    if len( dict_value ) == 0 : return None, None # if an empty dictionary is given, return None
    key_min = next( dict_value.__iter__( ) ) # retrieve the a key and value from the dictionary as an initial key-value pair with 'max' value
    value_min = dict_value[ key_min ]
    for key in dict_value :
        value = dict_value[ key ]
        if value_min > value :
            key_min = key
            value_min = value
    return key_min, value_min


# In[ ]:


def DICTIONARY_Merge_Min( * l_dict ) :
    ''' 
    # 2020-12-06 19:33:59 
    merge dictionaries by putting minimum value for each key '''
    dict_merged = l_dict[ 0 ]
    for dict_to_be_merged in l_dict[ 1 : ] :
        for key in dict_to_be_merged :
            dict_merged[ key ] = min( dict_merged[ key ], dict_to_be_merged[ key ] ) if key in dict_merged else dict_to_be_merged[ key ]
    return dict_merged

def DICTIONARY_Merge_Max( * l_dict ) :
    '''
    # 2020-12-06 19:33:59 
    merge dictionaries by putting maximum value for each key '''
    dict_merged = l_dict[ 0 ]
    for dict_to_be_merged in l_dict[ 1 : ] :
        for key in dict_to_be_merged :
            dict_merged[ key ] = max( dict_merged[ key ], dict_to_be_merged[ key ] ) if key in dict_merged else dict_to_be_merged[ key ]
    return dict_merged


# In[ ]:


def DICTIONARY_Filter( dict_data, above = None, below = None, include = True ) :
    '''  # 2020-12-08 22:36:58 
    filter the key-value pairs in the given dictionary based on the value with the given threshold. 
    For example, when above = 10, below = 20, and include = True, return key-value pairs with 10 < value < 20
    '''
    if above is None and below is None : 
        return dict_data # if no threshold is given, return the given dictionary as-it-is
    dict_data_filtered = dict( )
    for key in dict_data :
        value = dict_data[ key ]
        
        if above is not None and below is not None :
            if value < below or above < value :
                if include :
                    dict_data_filtered[ key ] = value
            elif not include :
                dict_data_filtered[ key ] = value
        elif above is not None and below is None :
            if above < value :
                if include :
                    dict_data_filtered[ key ] = value
            elif not include :
                dict_data_filtered[ key ] = value
        elif above is None and below is not None :
            if value < below :
                if include :
                    dict_data_filtered[ key ] = value
            elif not include :
                dict_data_filtered[ key ] = value
    return dict_data_filtered


# In[ ]:


def COUNTER( l_values, dict_counter = None, ignore_float = True ) : # 2020-07-29 23:49:51 
    ''' Count values in l_values and return a dictionary containing count values. if 'dict_counter' is given, countinue counting by using the 'dict_counter'. if 'ignore_float' is True, ignore float values, including np.nan '''
    if dict_counter is None : dict_counter = dict( )
    if ignore_float : # if 'ignore_float' is True, ignore float values, including np.nan
        for value in l_values :
            if isinstance( value, float ) : continue # ignore float values
            if value in dict_counter : dict_counter[ value ] += 1
            else : dict_counter[ value ] = 1
    else : # faster counting by not checking type of value
        for value in l_values :
            if value in dict_counter : dict_counter[ value ] += 1
            else : dict_counter[ value ] = 1
    return dict_counter


# In[ ]:




# ### Utility functions

# In[ ]:


def UTIL_Assign_integer_rank_using_hashing( arr ) :
    ''' Assign integer to each unique entries in a given array and return a dictionary mapping the entry to the integer number '''
    return pd.Series( dict( ( entry, entry.__hash__( ) ) for entry in set( arr ) ) ).argsort( ).to_dict( )
UTILITY_Assign_integer_rank_using_hashing = UTIL_Assign_integer_rank_using_hashing


# In[ ]:


def UTIL_PD_Display_a_row( df, int_index = 0 ) :
    ''' Return a numpy array of a row along with the column labels ''' 
    return np.vstack( ( df.columns.values, df.iloc[ int_index ].values ) ).T


# In[ ]:


def UTIL_PANDAS_Search_Columns( df, query ) :
    return list( col for col in df.columns.values if query in col )


# ### Functions for Glob.glob module

# In[ ]:


def GLOB_Retrive_Strings_in_Wildcards( str_glob, l_dir_match = None, return_dataframe = True, retrieve_file_size = False, retrieve_last_modified_time = False, time_offset_in_seconds = 3600 * 9 ) : # 2020-11-16 18:20:52 
    """ retrieve strings in '*' wildcards in list of matched directories for the given string containing '*' wildcards. return strings in wildcards as a nested lists. Consecutive wildcards should not be used ('**' should not be used in the given string)
    'retrieve_file_size': if 'return_dataframe' is True, return file sizes in bytes by using os.stat( dir_match ).st_size
    'retrieve_last_modified_time': return the last modified time with pandas datetime datatype
    'time_offset_in_seconds': offset in seconds to Coordinated Universal Time (UTC) """
    l_dir_match = glob.glob( str_glob ) if l_dir_match is None else l_dir_match # retrive matched directories using glob.glob if 'l_dir_match' is not given
    l_intervening_str = str_glob.split( '*' ) # retrive intervening strings in a glob string containing '*' wildcards 
    l_l_str_in_wildcard = list( )
    for dir_match in l_dir_match : # retrive strings in wildcards for each matched directory
        dir_match_subset = dir_match.split( l_intervening_str[ 0 ], 1 )[ 1 ]
        l_str_in_wildcard = list( )
        for intervening_str in l_intervening_str[ 1 : ] : 
            if len( intervening_str ) > 0 : str_in_wildcard, dir_match_subset = dir_match_subset.split( intervening_str, 1 )
            else : str_in_wildcard, dir_match_subset = dir_match_subset, '' # for the wildcard at the end of the given string, put remaining string into 'str_in_wildcard' and empties 'dir_match_subset'
            l_str_in_wildcard.append( str_in_wildcard )
        l_l_str_in_wildcard.append( l_str_in_wildcard )
    if return_dataframe : # return dataframe containing strings in wildcards and matched directory
        df = pd.DataFrame( l_l_str_in_wildcard, columns = list( 'wildcard_' + str( index ) for index in range( str_glob.count( '*' ) ) ) )
        df[ 'dir' ] = l_dir_match
        if retrieve_file_size : df[ 'size_in_bytes' ] = list( os.stat( dir_match ).st_size for dir_match in l_dir_match )
        if retrieve_last_modified_time : 
            df[ 'time_last_modified' ] = list( datetime.datetime.utcfromtimestamp( os.path.getmtime( dir_file ) + time_offset_in_seconds ).strftime( '%Y-%m-%d %H:%M:%S' ) for dir_file in df.dir.values )
            df.time_last_modified = pd.to_datetime( df.time_last_modified ) # convert to datetime datatype
        return df
    else : return l_l_str_in_wildcard


# ### Basic functions 

# ##### Writing and Reading Pickle files

# In[ ]:


def PICKLE_Write( dir_file, data_object ) :
    ''' write binary pickle file of a given data_object '''
    with open( dir_file, 'wb' ) as handle :
        pickle.dump( data_object, handle, protocol = pickle.HIGHEST_PROTOCOL )
def PICKLE_Read( dir_file ) :
    ''' write binary pickle file of a given data_object '''
    with open( dir_file, 'rb' ) as handle : 
        data_object = pickle.load( handle ) 
    return data_object


# ##### String functions

# In[ ]:


# Moved to STR module
"""def STR_Replace_str_index( text, index, replacement ):
    return '{}{}{}'.format( text[ : index ], replacement , text[ index + 1 : ] )
def STR_Insert_new_line_fasta( fasta_seq ) : # insert new line characters in the sequences every 60 characters
    return '\n'.join( list( fasta_seq[ 60 * index_line : 60 * ( index_line + 1 ) ] for index_line in np.arange( int( ( len( fasta_seq ) - 1 ) / 60 ) + 1 ) ) ) # Since there should not be additional no new additional character at the end for the sequences with 60 amino acid, ( len( seq ) - 1 ) / 60 will be used"""


# ##### Functions related with data types

# In[ ]:


def TYPE_Convert_NP_Array( data, dtype = None ) :
    '''Default dtype = Float'''
    if dtype is None :
        dtype = float
    if type( data ) is pd.Series :
        data = data.values.astype( dtype )
    elif type( data ) is list :
        data = np.array( data, dtype = dtype )
    elif type( data ) is set :
        data = np.array( list( data ), dtype = dtype )
    elif type( data ) is not np.ndarray :
        print( 'ERROR: Invalid data type' )
        return -1
    return data


# ##### Functions for mapping

# In[ ]:


# Moved to MAP module
"""def Map_a2b( a ) :
    if a in dict_a2b :
        return dict_a2b[ a ] 
    else :
        return np.nan 
def Remove_version_info( ID ) :
    return ID.split( '.' )[ 0 ]"""


# In[ ]:


def ENSEMBL_remove_version_info( list_ensembl_ID ) :
    ''' remove version info after '.' of ensembl_ID from a given list of ensembl_IDs '''
    return list( ensembl_ID.split( '.' )[ 0 ] if '.' in ensembl_ID else ensembl_ID for ensembl_ID in list_ensembl_ID )


# In[ ]:


def DICT_Subset_dictionary( dict_full, l_keys ) :
    ''' Subset dictionary with the given list of keys ("l_keys") '''
    return { a_key : dict_full[ a_key ] for a_key in dict_full.keys( ) & set( l_keys ) }

"""def DICT_Subset_dictionary( dict_full, l_keys ) : # 20190921
    ''' Subset dictionary with the given list of keys ("l_keys") '''
    dict_subset = dict( )
    for a_key in l_keys :
        if a_key in dict_full :
            dict_subset[ a_key ] = dict_full[ a_key ]
    return dict_subset"""


# In[ ]:


def LIST_Deduplicate( a_list ) :
    ''' Remove duplicates in a list while preserving the order '''
    seen = set( )
    seen_add = seen.add
    return list( entry for entry in a_list if not ( entry in seen or seen_add( entry ) ) )


# In[ ]:


def LIST_intersection_with_set( a_list, a_set ) :
    ''' Return intersection between a list and a set (the set can be list or np.array type), while preserving the order '''
    list_intersection = list( an_element for an_element in a_list if an_element in a_set )
    if type( a_list ) is np.ndarray : # if received array is numpy array, also return numpy array
        return np.array( list_intersection, dtype = object )
    else :
        return list_intersection

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
    
    
# In[ ]:


def ARRAY_sort_one_list_based_on_another( list_A, list_B ) :
    '''  Sort entries in list_A based on values in list_B, and return np.array of sorted list_A and list_B. The two lists have to be the same length  '''
    list_A, list_B = np.array( list_A ), np.array( list_B )
    list_B_argsort = list_B.argsort( )
    return list_A[ list_B_argsort ], list_B[ list_B_argsort ]


# In[ ]:


def To_python_compatible_str( a_string, dict_replacement = { '%' : '_Percent_', '+' : '_Plus_', '-' : '_Minus_', '&' : '_and_', '=' : '_Equal_' }  ) :
    ''' convert a string into python-compatible string. '''
    l_incompatible_character = [ ' ', '?', '(', ')', '&', '%', '/', ',', ':', '.', '-', '+', '[', ']', '#', '=', '\n', '"', '\\', '|', '?', '*' ]
    for incompatible_character in l_incompatible_character : 
        a_string = a_string.replace( incompatible_character, dict_replacement.get( incompatible_character, '_' ) )
    return a_string


# In[ ]:


def To_window_path_compatible_str( a_string ) :
    '''
        replace following characters to '_' so that a given string will be compatible for Window file system :
    : (colon)    " (double quote)    / (forward slash)    \ (backslash)    | (vertical bar or pipe)    ? (question mark)    * (asterisk)
        Also, replace new line character into '_'
    '''
    return a_string.replace( '\n', '_' ).replace( ':', '_' ).replace( '"', '_' ).replace( '/', '_' ).replace( '\\', '_' ).replace( '|', '_' ).replace( '?', '_' ).replace( '*', '_' )


# In[ ]:


def Search_list_of_strings( list_of_strings, query = 'cancer', return_mask_matched = False, return_location_matched = False ) :
    ''' search list of strings to find strings that contains query string and return the result as a list. if 'return_mask_matched' is True, 
    return list_mask for locations of matched entries (return np.array( search_result, dtype = object ), list_mask_matched)  '''
    search_result, list_mask_matched = list( ), list( )
    list_of_strings = list_of_strings.values if type( list_of_strings ) is pd.Series else list_of_strings # if type of list_of_strings is pd.Series, even though iterating through pandas series is just fine, to be safe and fast, convert pd.Series to a numpy array
    for string in list_of_strings : # search list of strings to find strings that contains query string and return the result as a list 
        if not isinstance( string, ( float, int ) ) and query in string :
            search_result.append( string )
            list_mask_matched.append( True )
        else :
            list_mask_matched.append( False )
    if return_mask_matched :
        return np.array( search_result, dtype = object ), np.array( list_mask_matched, dtype = bool )
    elif return_location_matched :
        return np.where( np.array( list_mask_matched, dtype = bool ) )[ 0 ]
    else :
        return search_result


# In[1]:


def Search_list_of_strings_with_multiple_query( l_str, * l_query, return_mask = False, return_position = False ) :
    ''' Search list of strings with multiple query. for negative query, add '-' in front of the query '''
    arr_mask_matched = np.ones( len( l_str ), dtype = bool )
    for query in l_query :
        bool_query_positive = True
        if query[ 0 ] == '-' :
            bool_query_positive, query = False, query[ 1: ]
        l_mask_matched_for_a_query = list( True if query in entry else False for entry in l_str ) if bool_query_positive else list( False if query in entry else True for entry in l_str )
        arr_mask_matched = arr_mask_matched & np.array( l_mask_matched_for_a_query, dtype = 'bool' )
    if return_position : return np.where( arr_mask_matched )[ 0 ]
    return arr_mask_matched if return_mask else np.array( l_str, dtype = object )[ arr_mask_matched ]


# In[ ]:


def Search_list_of_strings_Return_mask( data, query, is_negative_query = False ) :
    if is_negative_query :
        return np.array( list( False if query in entry else True for entry in data ), dtype = bool )
    else :
        return np.array( list( True if query in entry else False for entry in data ), dtype = bool )


# ##### Functions for Handling Sets

# In[1]:


def SET_A_B_Intersection( arr_A, arr_B, return_intersection = False, label_A = 'A', label_B = 'B' ) :
    set_A, set_B = set( arr_A ), set( arr_B )
    print( "{label_A} = {n_A}, {label_B} = {n_B}, {label_A} ∩ {label_B} = {n_intersection}".format( label_A = label_A, label_B = label_B, n_A = len( set_A ), n_B = len( set_B ), n_intersection = len( set_A.intersection( set_B ) ) ) )
    if return_intersection :
        return set_A.intersection( set_B )


# In[ ]:


def SET_A_B_union_intersec_diff( A, B ) :
    A, B = set( A ), set( B )    
    print( 'number of entries : {}, {}'.format( len( A ), len( B ) ) )
    print( 'union, intersection : {}, {}'.format( len( A.union( B ) ), len( A.intersection( B ) ) ) )
    print( 'difference of A and B : {}, {}'.format( len( A.difference( B ) ), len( B.difference( A ) ) ) )


# In[ ]:


def SET_union_list_sets( list_sets ) :
    ''' return a union (all unique elements) of sets in a given list of gene_sets  '''
    set_union = set( )
    for a_set in list_sets :
        set_union.update( a_set )
    return set_union


# In[ ]:


def LIST_Split( l = None, n_split = 0, return_slice = False ) : # 2020-07-30 01:04:31 
    """ split a list into 'n_split' number of chunks. if 'return_slice' is True, return slice() instances instead of actually spliting the given list-like object """
    if return_slice : return list( slice( index_split, None, n_split ) for index_split in np.arange( n_split ) ) # if 'return_slice' is True, return slice() instances instead of actually spliting the given list-like object
    else : return list( l[ index_split : : n_split ] for index_split in np.arange( n_split ) )


# In[ ]:


def pl( arr ):
    print( len ( arr ) ) 


# In[ ]:


def pls( arr ):
    print( len( set( arr ) ) )


# In[ ]:


def n_dup( arr ) :
    print( len ( arr ) - len( set( arr ) ) )


# In[ ]:


def print_tree( tree, print_elements = False, d = 0 ):
    '''
    recursive method to visualize a dictionary-based tree
    '''
    if (tree == None or len(tree) == 0):
        print ( "\t" * d, "-" )
    else:
        for key, val in tree.items():
            if (isinstance(val, dict)):
                print ( "\t" * d, key )
                print_tree( val, print_elements, d + 1 ) # recursion
            else:
                if type( val ) is set or type( val ) is list : # if leaf node is set or list type, print the number of items
                    if print_elements :
                        print ( "\t" * d, key )
                        for element in val :
                            print( "\t" * ( d + 1 ), element)
                    else :
                        print ( "\t" * d, key, '( ' + str( len( val ) ) + ' )' )
                else :
                    print ( "\t" * d, key, '( ' + val + ' )' )


# In[ ]:


def PHOSPHOPROTEOME_Gene_phosphosites_display( Gene, list_samples = None ) :
    ''' display phosphosites of a gene along with number of valid samples for each phosphosite as a DataFrame '''
    Gene_ID = Gene_2_Gene_ID( Gene ) # retrive Gene_ID
    if not Gene_ID in df_phosphoproteome.Gene_ID.values : # Check whether the given Gene_ID is valid
        print( 'Invalid gene' )
        return
    if list_samples is None :
        list_samples = dict_type_sample[ 'normal_or_tumor' ]
    # acquire list of phosphosites and list of numbers of valid samples
    Gene_data = df_phosphoproteome[ df_phosphoproteome.Gene_ID == Gene_ID ][ list_samples ] # data of list_samples of df_phosphoproteome for a given Gene 
    indices = Gene_data.index.values.tolist( )
    list_num_samples = list( len( list_samples ) - np.sum( np.isnan( Gene_data ) , axis = 1 ).values )
    return pd.DataFrame( dict( indices = indices, number_of_samples = list_num_samples ) ).set_index( 'indices' ) # Create DataFrame and return


# In[ ]:


def Gene_2_Gene_ID( Gene ) :
    '''    check given Gene_Symbol or Gene_ID is valid and if valid, convert it into Gene_ID and return Gene_ID    '''
    if type( Gene ) == str :
        Gene = Gene.upper( ) 
        if Gene in dict_Symbol_ID_simple :
            return dict_Symbol_ID_simple[ Gene ]
        elif Gene.upper( ) in dict_Symbol_ID_simple : # try again after converting lower case characters into upper clase characters
            return dict_Symbol_ID_simple[ Gene.upper( ) ]
        else :
#             print( 'Gene_Symbol does not exist' )
            return -1
    elif Gene in dict_ID_Symbol_simple :
        return Gene
    else :
#         print( 'Gene_ID does not exist' )
        return -1
    return -1


# In[ ]:


def List_Gene__2__List_Gene_ID( Genes, return_mask_valid_genes = False ) :
    ''' convert Gene to Gene_ID using the method 'Gene_2_Gene_ID' return valid Gene_IDs and mask of valid genes in the original list of genes if 'return_mask_valid_genes' is set to True  '''
    Gene_IDs = np.array( list( Gene_2_Gene_ID( Gene ) for Gene in Genes ) )
    mask_valid_genes = Gene_IDs != -1
    Gene_IDs = Gene_IDs[ mask_valid_genes ]
    if return_mask_valid_genes :
        return Gene_IDs, mask_valid_genes
    else :
        return Gene_IDs


# In[ ]:


def GET_Gene_ID_and_Symbol_from_input_gene( Gene ) :
    '''    check given Gene_Symbol or Gene_ID is valid and if valid, convert it into Gene_ID and Gene_Symbol and return Gene_ID, Sumbol '''
    if type( Gene ) == str :
        if Gene in dict_Symbol_ID_simple :
            return dict_Symbol_ID_simple[ Gene ], Gene
        elif Gene.upper( ) in dict_Symbol_ID_simple :
            Gene_Symbol = Gene.upper( )
            return dict_Symbol_ID_simple[ Gene_Symbol ], Gene_Symbol
        else :
            print( 'Gene_Symbol does not exist' )
            return -1, -1
    elif Gene in dict_ID_Symbol_simple :
        return Gene, dict_ID_Symbol_simple[ Gene ]
    else :
        print( 'Gene_ID does not exist' )
        return -1, -1
    return -1, -1


# In[ ]:


def List_Gene_ID__2__List_Gene_Symbol( list_Gene_ID, add_gene_name = False, return_sorted = False, limit_gene_name_len = 35, return_mask_mapped = False ) :
    '''    convert List of Gene_IDs into List of Gene_Symbol, or optionally, Geme_Symbol + Gene_Name
    add_gene_name = False : if True, add Gene_Name in addition to Gene_Symbol. 
    If 'return_mask_mapped' is True, List_Symbol, List_mask_mapped     '''
    List_mask_mapped = list( ) # a mask that indicate an entry has been mapped
    List_Symbol = [ ]
    for Gene_ID in list_Gene_ID : # for each Gene_ID
        if Gene_ID not in dict_ID_Symbol_simple :
            List_mask_mapped.append( False )
        else :
            List_mask_mapped.append( True )
            if add_gene_name : # convert to Geme_Symbol + Gene_Name
                gene_name = dict_ID_Symbol[ Gene_ID ][ 1 ] # retrive gene name
                if len( gene_name ) > limit_gene_name_len : # if length of gene_name exceed limit_gene_name_len, cut and add '..'
                    gene_name = gene_name[ : limit_gene_name_len - 2 ] + '..'
                Gene_label = gene_name + ' (' + dict_ID_Symbol[ Gene_ID ][ 2 ] + ')' # "Gene_Name (Geme_Symbol)"
                List_Symbol.append( Gene_label )
            else : # just convert to Geme_Symbol
                List_Symbol.append( dict_ID_Symbol_simple[ Gene_ID ] )
    if return_sorted :
        return sorted( List_Symbol ) # return list of symbols
    elif return_mask_mapped :
        return np.array( List_Symbol, dtype = object ), np.array( List_mask_mapped, dtype = bool ) # return list of symbols and mask of mapped symbols
    else :
        return List_Symbol


# In[ ]:


def List_Gene_Symbol__2__List_Gene_IDs( list_Gene_Symbols ) :
    '''  convert List of Gene_Symbol into List of Gene_IDs for annotation purposes  '''
    List_ID = [ ]
    for Gene_Symbol in list_Gene_Symbols : # for each Gene_Symbol
        if Gene_Symbol in dict_Symbol_ID_simple : # when Gene_Symbol is valid
            List_ID.append( dict_Symbol_ID_simple[ Gene_Symbol ] )
        else :
            print( Gene_Symbol, 'not exist in the annotation dictionary' ) # if Symbol is invalid, report the symbol
    return List_ID # return list of symbols


# In[ ]:


def List_Gene_Symbol__2__List_Gene_IDs( list_Gene_Symbols, return_mask_mapped = False ) :
    '''    Convert List of Gene_Symbol into List of Gene_IDs, and return List_Symbols.   If 'return_mask_mapped' is True, return List_Symbols, List_mask_mapped     '''
    List_mask_mapped = list( ) # a mask that indicate an entry has been mapped
    List_IDs = [ ]
    for Gene_Symbol in list_Gene_Symbols : # for each Gene_ID
        if Gene_Symbol not in dict_Symbol_ID_simple :
            List_mask_mapped.append( False )
        else :
            List_mask_mapped.append( True )
            List_IDs.append( dict_Symbol_ID_simple[ Gene_Symbol ] )
    if return_mask_mapped :
        return np.array( List_IDs, dtype = object ), np.array( List_mask_mapped, dtype = bool ) # return list of symbols and mask of mapped symbols
    else :
        return List_IDs


# In[ ]:


def List_Gene_ID_Phosphosite__2__List_Gene_IDs( list_Gene_ID_Phosphosites ) :
    '''    Convert List of Gene_ID_Phosphosite of df_phosphoproteome into List of Gene_IDs, and return List_Gene_IDs   '''
    return list( float( Gene_ID_Phosphosite.split( '|' )[ 0 ] ) for Gene_ID_Phosphosite in list_Gene_ID_Phosphosites )


# In[ ]:


def SEARCH_gene_symbol_name( df = None, query_Symbol = None, query_Name = None, negative_query = None ) :
    ''' Using df (a subset of HGNC) or full HGNC df by default, search a given query and negative query, and return a DataFrame with genes that 
    have matching Gene_Symbol or Gene_Name to a given query. Give previous search result through an argument 'df' to refine search  '''
    if df is None : # set default df, which is a full HGNC table
        df = df_ID_Symbol
    if query_Symbol is not None :
        search_list, query = df.Approved_Symbol.values, query_Symbol
    elif query_Name is not None :
        search_list, query = df.Approved_Name.values, query_Name
    else :
        return - 1
    search_result, list_mask_matched = Search_list_of_strings( search_list, query = query, return_mask_matched = True )
    if negative_query is not None : # if negative query has been given, remove entries that were matched with the negative query
        search_result_neg_query, list_mask_matched_neg_query = Search_list_of_strings( search_list, query = negative_query, return_mask_matched = True )
        df_search_result = df.loc[ list_mask_matched & ( ~ list_mask_matched_neg_query ) ]
    else :
        df_search_result = df.loc[ list_mask_matched ]
    return df_search_result


# In[ ]:


def GET_list_Gene_IDs_Phosphosites__from__list_Gene_IDs( Gene_IDs, print_messages = False ) :
    '''     Convert a list of Entrez Gene_IDs into a list of df_phosphoproteome_indices for the use in df_phosphoproteome     '''
    list_phosphoproteome_indices = list( )
    for Gene_ID in Gene_IDs :
        if Gene_ID in dict_Gene_ID__list_phosphoproteome_indices :
            list_phosphoproteome_indices.extend( dict_Gene_ID__list_phosphoproteome_indices[ Gene_ID ] )
        elif print_messages :
            print( Gene_ID, 'not exist in df_phosphoproteome')
    return list_phosphoproteome_indices


# In[ ]:


def non_nan_percentage ( arr ) :
    '''
    print a percentage of non-NaN samples of an entry
    '''
    arr = np.array( arr, dtype = float ) # convert to float numpy array
    length_arr = len( arr ) # the length of array
    num_non_nan = np.sum( ~ np.isnan( arr ) ) # the number of non-nan elements in the array
    print( 'length :', length_arr, ', non-NaN percantage :', round( num_non_nan / length_arr * 100, 2 ), '%' )


# In[ ]:


def non_nan_percentages ( arr ) :
    '''
    return percentages of non-NaN samples of multiple entries as a numpy array
    '''
    arr = np.array( arr, dtype = float ) # convert to float numpy array
    length_arr = len( arr[ 0 ] ) # the length of array
    num_non_nan = np.sum( ~ np.isnan( arr ), axis = 1 ) # the number of non-nan elements in the array
    return( num_non_nan / length_arr )


# In[ ]:


def common_non_NaN_values( data_1, data_2 ) :
    '''    return data pairs with only non_NaN values    '''
    data_1 = np.array( data_1, dtype = float ) # convert array into numpy array
    data_2 = np.array( data_2, dtype = float )
    both_non_NaN = ( ~ np.isnan( data_1 ) ) & ( ~ np.isnan( data_2 ) )
    return data_1[ both_non_NaN ], data_2[ both_non_NaN ] # return data where only pairs with non NaN values are retained 


# In[ ]:


def get_data__tumor_grades( gene_id, df = None ) :
    '''
    return dictionary of tumor_grade-wise proteome data of a given gene 
    df = df_proteome : dataframe from which data will be retrived
    '''
    categories = [ 'G1', 'G2', 'G3', 'G4' ] # set categories
    # retrive data of each tumor_grade 
    G1_data = df[ list( dict_tumor_grade__sample[ categories[ 0 ] ] ) ].loc[ gene_id ].values
    G2_data = df[ list( dict_tumor_grade__sample[ categories[ 1 ] ] ) ].loc[ gene_id ].values
    G3_data = df[ list( dict_tumor_grade__sample[ categories[ 2 ] ] ) ].loc[ gene_id ].values
    G4_data = df[ list( dict_tumor_grade__sample[ categories[ 3 ] ] ) ].loc[ gene_id ].values
    return dict( G1 = G1_data, G2 = G2_data, G3 = G3_data, G4 = G4_data ) # return dictionary of data of genes 


# In[ ]:


def FIND_Gene_Sets__with_Gene_or_Gene_Set_Name( dict_Gene_Sets, Gene = None, Gene_Set_Name = None) :
    '''  Find and return a list of Gene_Sets in a dictionary of Gene_Sets that include a given Gene  '''
    list_matched_Gene_Sets = list( )
    if Gene is not None :
        Gene_ID = Gene_2_Gene_ID( Gene ) # convert Gene (Symbol or Gene_ID) to Gene_ID
        if Gene_ID == -1 :
            return -1
        for Gene_Set in dict_Gene_Sets : # for each Gene_Set in a dict_Gene_Sets, check whether a given Gene_ID exist in the Gene_Set
            if Gene_ID in dict_Gene_Sets[ Gene_Set ] :
                list_matched_Gene_Sets.append( Gene_Set )
    elif Gene_Set_Name is not None :
        for Gene_Set in dict_Gene_Sets : # for each Gene_Set in a dict_Gene_Sets, check whether a given Gene_ID exist in the Gene_Set
            if Gene_Set_Name in Gene_Set :
                list_matched_Gene_Sets.append( Gene_Set )       
    else :
        print( 'invalid query' )
        return -1
    return list_matched_Gene_Sets # return a list of Gene_Sets with a given Gene_IDs


# In[ ]:


def GET_series_sample_data_from_series_patient_data( series_patient_data ) :
    ''' Convert Series where index is patient_id to Series with sample_id '''
    dict_patient_data = series_patient_data.to_dict()
    dict_sample_data = dict( )
    for Sample_ID, Patient_ID in df_samples.loc[ All_samples ].ParticipantID.to_dict( ).items( ) :
        dict_sample_data[ Sample_ID ] = dict_patient_data[ Patient_ID ]
    series_sample_data = pd.Series( dict_sample_data )
    return series_sample_data.rename( series_patient_data.name )


# #### Functions for handling intervals

# In[ ]:


def INTERVAL_Overlap( interval_1, interval_2, flag_sort_to_retrieve_start_and_end = False, flag_0_based_coordinate_system = True ) : # 2020-08-06 20:44:47 
    ''' Fast, basic function for retrieving overlap between two intervals.
    return number of overlapped length between two intervals (each interval is a tuple or list containing start and end position). 
    'flag_0_based_coordinate_system': if interval contains float numbers, set 'flag_0_based_coordinate_system' to True.
    'flag_sort_to_retrieve_start_and_end': if interval is always (start, end), set 'flag_sort_to_retrieve_start_and_end' to False to improve performance. (about 200ns faster)  '''
    if flag_sort_to_retrieve_start_and_end :
        start_1, end_1 = sorted( interval_1 )
        start_2, end_2 = sorted( interval_2 )
    else :
        start_1, end_1 = interval_1
        start_2, end_2 = interval_2
    if not flag_0_based_coordinate_system : start_1, start_2 = start_1 - 1, start_2 - 1
    if ( end_1 <= start_2 ) | ( end_2 <= start_1 ) : return 0
    else : return min( end_1, end_2 ) - max( start_1, start_2 )


# #### Functions for handling outliers

# In[ ]:


def OUTLIERS_GET_mask_for_outliers( arr, n_std_for_outliers, outlier_percentile_std = 5 ) :
    ''' for each row of a given numpy array, identify the outliers (outside n_std_for_outliers * std from mean) and return a boolean array 
        that indicates locations of outliers. if 'outlier_percentile_std' is not zero, calculate and mean and std excluding top and bottom 'outlier_percentile_std' percentile for defining outliers  '''
    if outlier_percentile_std == 0 :
        arr_mean, arr_std = arr.mean( axis = 1 ), arr.std( axis = 1 ) # retrive mean and std of outliers
    else :
        arr_masked = np.ma.masked_array( data = arr, mask = ( arr.T < np.percentile( arr, outlier_percentile_std, axis = 1 ) ).T | ( arr.T > np.percentile( arr, 100 - outlier_percentile_std, axis = 1 ) ).T )
        arr_mean, arr_std = arr_masked.mean( axis = 1 ).data, arr_masked.std( axis = 1 ).data # retrive mean and std of outliers
    lower_limit, upper_limit = arr_mean - n_std_for_outliers * arr_std, arr_mean + n_std_for_outliers * arr_std # set limits for outliers
    return ( ( arr.T < lower_limit ) | ( upper_limit < arr.T ) ).T # return mask indicates location of outliers


# In[ ]:


def GET_non_NaN_without_outliers( data_1, data_2, n_std_for_outliers = 3 ) :
    ''' remove non_NaN pairs and pairs with outlier data values from a pair of arrays, data_1 and data_2, and return processed pair of arrays '''
    data_1, data_2 = common_non_NaN_values( data_1, data_2 ) # build remove NaN values
    mask_data_1, mask_data_2 = OUTLIERS_GET_mask_for_outliers( np.array( [ data_1, data_2 ] ), n_std_for_outliers = n_std_for_outliers ) # get mask for outliers for each gene
    mask_non_outliers = ~ ( mask_data_1 | mask_data_2 ) # build a mask for pairs of data that do not contain outliers
    return data_1[ mask_non_outliers ], data_2[ mask_non_outliers ] # return data without outliers


# ### Functions using Masked arrays

# In[ ]:


def MA_Mask_Zeros( df ) :
    data = df.replace( 0, np.nan ).values
    data = np.ma.masked_array( data, np.isnan( data ) )
    return data


# #### Functions for Series and Dataframes

# ##### Basic functions for Series

# In[ ]:


def S_percentile_bottom_top( s, bottom_percentile = 5, top_percentile = 95, return_label = True ) :
    ''' return subset of series of thresholded by the given percentiles. return only labels of subsets of top and bottom 5% by default '''
    s_bottom, s_top = s[ s < np.percentile( s, bottom_percentile ) ], s[ s > np.percentile( s, top_percentile ) ]
    return ( s_bottom.index.values, s_top.index.values ) if return_label else ( s_bottom, s_top )


# ##### Basic functions for DataFrames

# In[ ]:


def OPTION_PANDAS_Display_n_rows( n_rows = None ) :
    ''' Set number of rows to display. set 'n_rows' to None to reset the option,  '''
    if n_rows is not None :
        pd.set_option( 'display.max_rows', n_rows )
    else :
        pd.reset_option( 'display.max_rows' )


# In[ ]:


def PANDAS_Reindex( df ) :
    ''' set default index of a given dataframe of series '''
    return df.rename( index = dict( ( index_before, index_after ) for index_before, index_after in zip( df.index.values, np.arange( len( df ) ) ) ) )


# In[ ]:


def PANDAS_filter_by_sum_of_row( df, thres = None ) :
    arr_row_sum = df.values.sum( axis = 1 )
    if thres is None :
        thres = np.percentile( arr_row_sum, 95 ) # set 95% percentil as a threshold
    return df[ arr_row_sum > thres ]


# In[ ]:


def PANDAS_DISPLAY_list_df_loc( list_df, indices, return_sliced_df = False ) :
    ''' display list of df with a given index '''
    for df in list_df :
        display( df.loc[ indices ] )
    if return_sliced_df :
        return list( df.loc[ indices ] for df in list_df )


# In[ ]:


def PANDAS_DATAFRAME_number_unique_entries( df ) :
    ''' return a series with column labels containing number of unique entries excluding NaN value for each column label  '''
    dict_column_n_unique = dict( )
    for column in df.columns.values :
        dict_column_n_unique[ column ] = len( df[ column ].dropna( ).unique( ) )
    return pd.Series( dict_column_n_unique )


# In[ ]:


def DF_Build_Index_Using_Dictionary( df, l_col_for_index ) : # 2020-08-06 17:12:59 
    ''' return a dictionary with key = index or multi-index for a given 'l_col_for_index' and value = list of integer index for df.values (not df.index.values)
    Using Python dictionary and numpy array can be upto ~2000 times faster than using pandas.DataFrame.loc[ ] method for accessing multi-indexed rows. '''
    dict_index = dict( )
    if isinstance( l_col_for_index, str ) : # when only single col_name was given for index
        for int_index, index in enumerate( df[ l_col_for_index ].values ) :
            if index in dict_index : dict_index[ index ].append( int_index )
            else : dict_index[ index ] = [ int_index ]
    else : # when a list of col_names was given for index
        for int_index, arr_index in enumerate( df[ l_col_for_index ].values ) :
            t_index = tuple( arr_index )
            if t_index in dict_index : dict_index[ t_index ].append( int_index )
            else : dict_index[ t_index ] = [ int_index ]
    return dict_index 


# In[ ]:


def DF_Count_and_Drop_Duplicates( df, l_col_for_identifying_duplicates, l_col_for_sorting_and_not_for_identifying_duplicates = None, dict_kw_sort_values = { 'ascending' : False }, col_name_for_duplicate_counts = 'duplicate_counts', inplace = False, dict_kw_drop_duplicates = { 'keep' : 'first' } ) : # 2020-08-08 19:58:43  # pandas >= 1.1
    ''' drop duplicates while counting duplicates. if a column for counting duplicates already exist in the given dataframe, sum values in the columns while dropping duplicates instead of counting duplicates.
    'l_col_for_sorting_and_not_for_identifying_duplicates': if 'l_col_for_sorting_and_not_for_identifying_duplicates' is given, sort rows using columns given in the order in 'l_col_for_sorting_and_not_for_identifying_duplicates' before dropping duplicates '''
    if not inplace : df = deepcopy( df )
    if l_col_for_sorting_and_not_for_identifying_duplicates is not None : df.sort_values( l_col_for_identifying_duplicates + l_col_for_sorting_and_not_for_identifying_duplicates, inplace = True, ** dict_kw_sort_values ) # 'l_col_for_sorting_and_not_for_identifying_duplicates': if 'l_col_for_sorting_and_not_for_identifying_duplicates' is given, sort rows using columns given in the order in 'l_col_for_sorting_and_not_for_identifying_duplicates' before dropping duplicates
    if not isinstance( l_col_for_identifying_duplicates, ( list ) ) : l_col_for_identifying_duplicates = list( l_col_for_identifying_duplicates ) 
    set_col_for_identifying_duplicates = set( l_col_for_identifying_duplicates )
    l_col = list( col for col in df.columns.values if col != col_name_for_duplicate_counts ) # retrieve all column names except for a column for duplicate counts
    df_containing_col_not_for_identifying_duplicates = df.drop_duplicates( subset = l_col_for_identifying_duplicates, ** dict_kw_drop_duplicates ).set_index( l_col_for_identifying_duplicates ) # this dataframe can contains a column for duplicate counts
    if col_name_for_duplicate_counts in df.columns.values : df_containing_col_not_for_identifying_duplicates.drop( columns = col_name_for_duplicate_counts, inplace = True )
    else : df[ col_name_for_duplicate_counts ] = 1 # if a column for counting duplicates already exist in the given dataframe, sum values in the columns while dropping duplicates instead of counting duplicates.
    df = df[ l_col_for_identifying_duplicates + [ col_name_for_duplicate_counts ] ].groupby( l_col_for_identifying_duplicates, dropna = False ).sum( ).join( df_containing_col_not_for_identifying_duplicates ) # count duplicates by summing 1 or the previous duplicate counts given by the 'col_name_for_duplicate_counts' column
    df.reset_index( inplace = True )
    df = df[ l_col + [ col_name_for_duplicate_counts ] ] # rearange columns to match the order in the original dataframe
    return df


# In[ ]:


def DF_COLUMN_labels_replace_incompatible_characters( df, dict_replacement = { '%' : '_Percent_', '+' : '_Plus_', '-' : '_Minus_', '&' : '_and_', '=' : '_Equal_', '/' : '_by_' } ) :
    ''' remove incompatible characters in column labels of a given DataFrame so that columns are accesible without using [ ] bracket '''
    incompatible_characters = [ ' ', '?', '(', ')', '&', '%', '/', ',', ':', '.', '-', '+', '[', ']', '#', '=', '\n', '"', '\\', '|', '?', '*', "'", '−' ] # define incompatible characters
    dict_column_b4_column_after = dict( )  # an empty dictionary that will be used to rename columns
    columns = df.columns.values # retrive current columns
    for column in columns : # for each column
        column_after = str( column ).strip( ) # copy current column label            
        for incompatible_character in incompatible_characters : # for each incompatible charactor, if the charactor exists in column, replace it with '_' 
            if incompatible_character in column_after :
                column_after = column_after.replace( incompatible_character, dict_replacement.get( incompatible_character, '_' ) )
        if column_after[ 0 ].isdigit( ) : # if first character is a digit character, add '_' in front of the label 
            column_after = '_' + column_after
        dict_column_b4_column_after[ column ] = column_after # build a dictionary for renaming the columns
    return df.rename( columns = dict_column_b4_column_after ) # return dataFrame with changed column_labels
PANDAS_COLUMN_labels_replace_incompatible_characters = DF_COLUMN_labels_replace_incompatible_characters # define aliases for this function


# In[ ]:


def DF_Summary_CV( df, add_gene_annotation = True ) :
    ''' Summarize a given dataframe (preferably gene-expression matrix) '''
    df_summary = pd.DataFrame( dict( AVG = df.mean( axis = 1 ), STD = df.std( axis = 1 ) ) )
    df_summary[ 'CV' ] = df_summary.STD / df_summary.AVG
    df_summary = df_summary.sort_values( 'CV' )
    return PD_Add_gene_annotation( df_summary ) if add_gene_annotation else df_summary


# In[ ]:


def GET_MASK_of_intersection( arr, annotation ) :
    ''' return a boolean mask with shape of 'arr' in which True means the entry belongs to 'annotation'  '''
    list_mask = list( True if entry in annotation else False for entry in arr ) # create mask in list type
    return np.array( list_mask, dtype = bool ) # convert it to numpy boolean type array 


# In[ ]:


def GET_MASK_of_intersection( arr, annotation ) :
    ''' return a boolean mask with shape of 'arr' in which True means the entry belongs to 'annotation'  '''
    list_mask = list( True if entry in annotation else False for entry in arr ) # create mask in list type
    return np.array( list_mask, dtype = bool ) # convert it to numpy boolean type array 



def PD_Threshold( df, AND_operation = True, ** dict_thresholds ) :
    '''  Select rows of a given DataFrame or indices of Series based on a given threshold for each given column or the given series. 
    Add 'b' or 'B' at the end of column_label to select rows below the threshold, or add 'a' or 'A' to select rows above the threshold.
    If 'AND_operation' is true, filter generated from the given threshold will be combined with AND operation before filtering rows of a given dataframe  '''
    set_df_columns = set( df.columns.values ) if type( df ) is pd.DataFrame else set( [ '' ] )
    mask_filter = np.ones( len( df ), dtype = bool ) if AND_operation else np.zeros( len( df ), dtype = bool )
    for col_direction, threshold in dict_thresholds.items( ) :
        col, direction = col_direction[ : -1 ], col_direction[ -1 ]
        if col not in set_df_columns :
            print( "'{}' column_label does not exist in the given DataFrame".format( col ) )
            continue
        data = df[ col ].values if type( df ) is pd.DataFrame else df.values
        if direction.lower( ) == 'a' :
            current_mask = ( data > threshold )
        elif direction.lower( ) == 'b' :
            current_mask = ( data < threshold )
        else :
            print( "'{}' direction is not either 'a' or 'b' and thus invalid".format( direction ) )
            continue
        mask_filter = current_mask & mask_filter if AND_operation else current_mask | mask_filter
    return df[ mask_filter ]


# In[ ]:


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

# In[ ]:


def PD_Subset( df, subset = None, axis = 0, index = None, columns = None, preserve_order_in_df = True, subset_difference = False, return_numeric_column_only = False ) :
    ''' 
    # 2021-02-01 16:10:25 
    Subset a dataframe or series with entries given by 'subset' while preserving order in df (by not using set operations).
    if 'preserve_order_in_df' is True, preserving the order in df after subset. If False, preserve the order in a given subset . 
    If 'axis' is 0, subset indices, while 'axis' is 1, subset columns
    if 'subset_difference' is True, exclude entries in a given subset 
    'return_numeric_column_only' : if True, return columns with only numeric datatypes
    '''
    if index is not None : # subset df according to given index and columns before subsetting
        df = Subset( df, subset = index, axis = 0, preserve_order_in_df = preserve_order_in_df )
    if columns is not None :
        df = Subset( df, subset = columns, axis = 1, preserve_order_in_df = preserve_order_in_df )
    if subset is not None :
        if type( subset ) is pd.Series :
            subset = LIST_Deduplicate( subset.values )
        else :
            subset = LIST_Deduplicate( subset ) # convert list of values into set
        if type( df ) is pd.Series : # if pandas Series is given, only indices will be subsetted since columns do not exist
            axis = 0
        set_subset, set_index = set( subset ), set( df.index.values )
        if axis == 0 :
            if subset_difference : 
                df = df[ list( False if index in set_subset else True for index in df.index.values ) ] if preserve_order_in_df else df.loc[ list( entry for entry in subset if entry not in set_index ) ]
            else :
                df = df[ list( True if index in set_subset else False for index in df.index.values ) ] if preserve_order_in_df else df.loc[ list( entry for entry in subset if entry in set_index ) ]
        elif axis == 1 :
            set_columns = set( df.columns.values )
            if subset_difference : 
                df = df.loc[ :, list( False if column in set_subset else True for column in df.columns.values ) ] if preserve_order_in_df else df[ list( entry for entry in subset if entry not in set_columns ) ]
            else :
                df = df.loc[ :, list( True if column in set_subset else False for column in df.columns.values ) ] if preserve_order_in_df else df[ list( entry for entry in subset if entry in set_columns ) ]
    if return_numeric_column_only and len( df ) > 0 : # 'return_numeric_column_only' : if True, return columns with only numeric datatypes
        df = df[ list( key for key, value in df.iloc[ 0 ].to_dict( ).items( ) if isinstance( value, ( int, float ) ) ) ]
    return df


# In[ ]:


def PD_Search( df, query_AND_operation = True, is_negative_query = False, ignore_case = True, ** dict_search ) :
    ''' to search index, put 'index' as the name of column in the 'dict_search' argument. To search values in a series, put 'values' as the name of column in the 'dict_search' argument
    To retrive entries matched with all queries, set 'query_AND_operation' to True (perform AND operation during searching), while set 'query_AND_operation' to False to retrive all entries 
    containing at least one of the given queries. Search with multiple columns are always performed with AND operation (will be improved) '''
    for col, queries in dict_search.items( ) :
        if type( df ) is pd.DataFrame and col not in list( df.columns.values ) + [ 'index', 'columns', 'col' ] : # print out an error message when invalid column name was given for a dataframe
            print( "'{}' does not exist in columns of a given DataFrame".format( col ) )
            continue
        if col == 'index' : # retrive data to be searched
            data = df.index.values
        elif col == 'columns' or col == 'col' :
            data = df.columns.values
        elif col == 'values' :
            if type( df ) is not pd.Series : continue # if df is not Series, ignore this query
            data = df.values
        else : data = df[ col ].values
        if ignore_case : data = np.char.lower( data.astype( str ) ).astype( object ) # convert all strings to lower cases if 'ignore_case' is set to True.
        else : data = data.astype( str ).astype( object ) # convert data_types of data values to be searched to string
#         if isinstance( data[ 0 ], ( float, int, np.float64, np.int64 ) ) : # if data contained in the column (or values or indices) are numbers and not a string, ignore this column
#             continue
        if isinstance( queries, ( list, tuple, np.ndarray, set ) ) : # if the given queries is iterable
            len_entries = df.shape[ 1 ] if col == 'columns' or col == 'col' else df.shape[ 0 ]
            mask = np.ones( len_entries ).astype( bool ) if query_AND_operation else np.zeros( len_entries ).astype( bool ) # mask that will contain all entries matching a given list of data
            for query in queries :
                if ignore_case : query = query.lower( ) # convert query to lower cases if 'ignore_case' is set to True.
                mask_current_query = Search_list_of_strings_Return_mask( data, query, is_negative_query = is_negative_query )
                mask = mask & mask_current_query if query_AND_operation else mask | mask_current_query # merge mask for each data entry
        else :
            if ignore_case : query = queries.lower( ) # convert query to lower cases if 'ignore_case' is set to True.
            mask = Search_list_of_strings_Return_mask( data, query, is_negative_query = is_negative_query )
        df = df.iloc[ :, mask ] if col == 'columns' else df[ mask ]
    return df

# In[ ]:


def PD_Display( df, max_rows = 1000, max_columns = 35, max_colwidth = 200 ) :
    """ display more rows and columns of pandas objects. Also, display more column widths """
    with pd.option_context( 'display.max_rows', max_rows, 'display.max_columns', max_columns, 'display.max_colwidth', max_colwidth ) :
        display( df )


# In[ ]:


def PD_Column_Add_Suffix( df, suffix = '', inplace = False ) :
    """
    # 2021-02-01 16:37:41 
    add suffix to the columns
    """
    if len( suffix ) > 0 :
        if inplace : 
            df = deepcopy( df )
        df.columns = df.columns.values.astype( object ) + suffix
    return df


# In[ ]:

def PANDAS_Align_two( df_1, df_2, axis = 1 ) :
    ''' ** There is a original method "pandas.align" Use it if possible. **
    align two dataframe or series. by default, align by columns (axis = 1). if series is given as one of the inputs, align two objects by indices '''
    if type( df_1 ) is pd.Series or type( df_2 ) is pd.Series : # if series is given as one of the inputs, align two objects by indices
        axis = 0
    if axis == 1 :
        cols_shared = LIST_intersection_with_set( df_1.columns.values, df_2.columns.values )
        return df_1[ cols_shared ], df_2[ cols_shared ]
    else :
        indices_shared = LIST_intersection_with_set( df_1.index.values, df_2.index.values )
        return df_1.loc[ indices_shared ], df_2.loc[ indices_shared ]  


# In[ ]:


def PANDAS_Subset_relation_df( df, indices ) :
    ''' Return a subset of relation dataframe with Genes in the Gene_Set. A relation DataFrame is a n x n matrix where indices and columns in the same order. '''
    valid_indices = list( True if index in indices else False for index in df.index.values )
    return df.loc[ valid_indices, valid_indices ]


# In[ ]:


def PANDAS_relation_dataframe_make_symmetric( df, data_is_in_lower_tri = True, symmetry_relationship = 'same', diagonal_data = None ) :
    ''' For a given relation dataframe, copy data in lower triangle to upper triangle if 'data_is_in_lower_tri' is True (or copy data in upper triangle into lower triangle if 'data_is_in_lower_tri' 
    is False). 'symmetry_relationship' = 'same' or '=' ; 'opposite_sign' or '-' ; 'inverse' or '1/' . 'diagonal_data' : data_value that will be used to the fill diagonal data. If None is given,
    use the first diagonal value as 'diagonal_data' value ( df.iloc[ 0, 0 ] ). If 'diagonal_data' = 'retain', use the initial diagonal data '''
    n = len( df )
    if n < 1 :
        return -1
    indices_diagonal = np.diag_indices( n ) # define diagnoal indices 
    if diagonal_data is None :
        diagonal_data = df.iloc[ 0, 0 ]
    elif diagonal_data == 'retain' : # if 'diagonal_data' is set to 'retain', save diagonal values 
        arr_diagonal_data = df.values[ indices_diagonal ]
    if data_is_in_lower_tri :
        indices_data, indices_no_data = np.tril_indices( n ), np.triu_indices( n )
    else :
        indices_data, indices_no_data = np.triu_indices( n ), np.tril_indices( n )
    df = deepcopy( df )
    df.values[ indices_no_data ] = 0
    if symmetry_relationship == 'same' or symmetry_relationship == '=' :
        df = df + df.T
    elif symmetry_relationship == 'opposite_sign' or symmetry_relationship == '-' :
        df = df - df.T
    elif symmetry_relationship == 'inverse' or symmetry_relationship == '1/' :
        df_T = deepcopy( df.T )
        df_T.values[ indices_data ] = 1
        df_T = 1.0 / df_T
        df_T.values[ indices_data ] = 0
        df = df + df_T
    else :
        return -1
    if diagonal_data == 'retain' : # if 'diagonal_data' is set to 'retain', put a saved diagonal values back to relation_data matrix 
        df.values[ indices_diagonal ] = arr_diagonal_data
    else :
        df.values[ indices_diagonal ] = diagonal_data      
    return df


# ##### Functions for exploratory analysis on Pandas dataframs  

# In[ ]:


def PDA_Calculate_CV( df ) :
    ''' Calculate coefficient of variation (std / average) for values of each rows '''
    df = deepcopy( df )
    df = df[ df.sum( axis = 1 ) != 0 ] # remove rows with total coutns == 0 sort rows based on averages
    df = df.iloc[ df.sum( axis = 1 ).values.argsort( )[ : : -1 ] ] # sort rows based on total counts of each row
    dfcv = pd.DataFrame( dict( AVG = df.mean( axis = 1 ), STD = df.std( axis = 1 ) ) ) # calculate summary of df
    dfcv[ 'CV' ] = dfcv.STD.values / dfcv.AVG.values
    return dfcv


# ### Pandas functions for df containing genes and analysis results

# In[ ]:


def PD_Locate_gene( df, gene ) :
    ''' Access gene in as a geneid in the indices of columns of a given df '''
    geneid = Gene_2_Gene_ID( gene )
    if geneid == -1 : # check whether the given gene is valid
        return df.loc[ gene ]
    else :
        if isinstance( df.index.values[ 0 ], ( int, float, np.float64, np.int64 ) ) : # first check indices of a given df # if index contain a number (potentially geneid)
            if geneid in df.index.values :
                return df.loc[ geneid ]
        elif type( df ) is pd.DataFrame : # if the type of indices is other than a number, check whether it is dataframe, and check whether geneid is in the columns
            if geneid in df.columns.values :
                return df[ geneid ]
    return -1


# In[ ]:


def PD_Add_gene_annotation( df, in_place = False, geneid_to_genesymbol = False, type_gene_id = 'ensembl' ) :
    ''' For dataframe and series whose indices are Entrez Gene_IDs or df with multiindex with their first index as Entrez Gene_IDs, add Gene_Symbol and Gene_Name columns and return the dataframe. 
    If 'geneid_to_genesymbol' is True, map gene_id in its index to gene_symbols.
    Currently, Entrez Gene_IDs and Ensembl Gene_IDs are supported. '''
    if not in_place : # if gene annotation is added not in_place, copy dataframe before adding annotations
        df = deepcopy( df )
    if geneid_to_genesymbol :
        df_annotated = PD_Add_gene_annotation( df ).set_index( 'Approved_Symbol' ).drop( columns = 'Approved_Name' )
        return df_annotated.iloc[ :, 0 ] if type( df ) is pd.Series else df_annotated # return series if series was given
    Gene_IDs = df.index.values
    if type( Gene_IDs[ 0 ] ) is tuple : # if df are indexed with multiIndex, retrive only the first index from index tuple as a list of Gene_IDs
        Gene_IDs = np.array( list( multiIndex[ 0 ] for multiIndex in Gene_IDs ), dtype = object ) 
#     elif type( Gene_IDs[ 0 ] ) is str and '|' in Gene_IDs[ 0 ] :  # if df are indexed with phosphosite annotation, retrive only Gene_IDs from the index
#         Gene_IDs = np.array( list( Gene_ID.split( '|' )[ 0 ] for Gene_ID in df_res.index.values ), dtype = float )
    df_anno = df_anno_biomaRt_gene if ( ( type( Gene_IDs[ 0 ] ) is str and 'ENS' in Gene_IDs[ 0 ] ) and type_gene_id == 'ensembl' ) else df_ID_Symbol
    set_valid_gene_id = set( df_anno.index.values )
    mask_mapped = list( True if gene_id in set_valid_gene_id else False for gene_id in Gene_IDs )
    if np.sum( mask_mapped ) != len( df ) : # if there are some entries that are not mapped with current annotation data, retrive only valid entries  
        Gene_IDs = Gene_IDs[ mask_mapped ]
        df = df[ mask_mapped ]
    df_annotation = df_anno.loc[ Gene_IDs ] # retrive Gene_Name and Gene_Symbols of valid Gene_IDs
    arr_gene_names = df_annotation.Approved_Name.values
    arr_gene_symbols = df_annotation.Approved_Symbol.values
    if type( df ) is pd.DataFrame :
        if 'Approved_Symbol' in df.columns.values : # return the given dataframe if the given df already has annotation columns 
            return df
        df[ 'Approved_Name' ], df[ 'Approved_Symbol' ] = arr_gene_names, arr_gene_symbols
        return df[ [ 'Approved_Symbol', 'Approved_Name' ] + list( df.columns.values )[ : -2 ] ] # reorder the columns and return the dataframe
    elif type( df ) is pd.Series :
        data_name = df.name if df.name is not None else 'unknown_data'
        return pd.DataFrame( dict( Approved_Symbol = arr_gene_symbols, Approved_Name = arr_gene_names, data = df.values ), index = Gene_IDs ).rename( columns = dict( data = data_name ) )


# In[ ]:


"""def PANDAS_add_gene_annotation_to_df( df ) : # frozen 20190602
    ''' For dataframe whose indices are Entrez Gene_IDs or df with multiindex with their first index as Entrez Gene_IDs, add Gene_Symbol and Gene_Name columns and return the dataframe '''
    columns = list( df.columns.values )
    Gene_IDs = df.index.values
    if type( Gene_IDs[ 0 ] ) is tuple : # if df are indexed with multiIndex, retrive only the first index from index tuple as a list of Gene_IDs
        Gene_IDs = np.array( list( multiIndex[ 0 ] for multiIndex in Gene_IDs ), dtype = object ) 
        _, mask_mapped = List_Gene_ID__2__List_Gene_Symbol( Gene_IDs, return_mask_mapped = True )
        if np.sum( mask_mapped ) != len( df ) : # if there are some entries that are not mapped with current annotation data, retrive only valid entries  
            Gene_IDs = Gene_IDs[ mask_mapped ]
            df = df[ mask_mapped ]
        df_annotation = df_ID_Symbol.loc[ Gene_IDs ] # retrive Gene_Name and Gene_Symbols of valid Gene_IDs
        df[ 'Approved_Name' ] = df_annotation.Approved_Name.values
        df[ 'Approved_Symbol' ] = df_annotation.Approved_Symbol.values
        return df[ [ 'Approved_Symbol', 'Approved_Name' ] + columns ] # reorder the columns and return the dataframe
    else :
        list_Gene_IDs_valid = LIST_intersection_with_set( Gene_IDs, df_ID_Symbol.index.values ) # retrive list of valid Gene_IDs
        df = df.loc[ list_Gene_IDs_valid ] # retrive only valid Gene_IDs 
        return df.join( df_ID_Symbol[ [ 'Approved_Name', 'Approved_Symbol' ] ].loc[ list_Gene_IDs_valid ] )[ [ 'Approved_Symbol', 'Approved_Name' ] + columns ] # add Gene_Name and Gene_Symbols to df, and rearrange the columns so that annotation columns are located at left"""


# In[ ]:


def RESULT_filter( df_result, thres_abs_value = 0, label_value = None, thres_p_value = 1, thres_adj_p_value = 1, label_adj_p_value = 'adjusted_p_value', label_p_value = 'p_value', thres_n_gene = 0, label_n_genes = 'number_valid_genes', thres_abs_log2fc = 0, label_log2fc = 'Log2_Fold_Change', thres_abs_diff = 0, label_diff = 'Difference', thres_avg_all = 0, label_avg_all = 'average_all', log2fc_positive = None ) :
    ''' Filter out entries to retrive only significant results from (mostly) Log2FC analysis result dataframe. if 'log2fc_positive' is True, return only positive log2fc, and negative log2fc if False. '''
    if thres_n_gene > 0 : # if valid thres_n_gene is given, drop entries with n_genes below 'thres_n_gene'
        df_result = df_result[ df_result[ label_n_genes ] >= thres_n_gene ]
    if thres_p_value < 1 :
        df_result = df_result[ df_result[ label_p_value ] <= thres_p_value ] # return p-value filtered result if valid thres_p_value is given
    if thres_adj_p_value < 1 :
        df_result = df_result[ df_result[ label_adj_p_value ] <= thres_adj_p_value ] # return p-value filtered result if valid thres_p_value is given
    if thres_abs_log2fc > 0 :
        df_result = df_result[ np.abs( df_result[ label_log2fc ] ) >= thres_abs_log2fc ] # return Log2FC filtered result if valid 'thres_abs_log2fc' is given (entries with absolute value of Log2FC below thres_abs_log2fc is filtered out)
    if thres_abs_diff > 0 :
        df_result = df_result[ np.abs( df_result[ label_diff ] ) >= thres_abs_diff ] # return difference filtered result if valid 'thres_abs_diff' is given (entries with absolute value of Difference below thres_abs_diff is filtered out)
    if thres_avg_all > 0 :
        df_result = df_result[ df_result[ label_avg_all ] >= thres_avg_all ]
    if thres_abs_value > 0 :
        df_result = df_result[ np.abs( df_result[ label_value ] ) <= thres_abs_value ]
    if log2fc_positive is not None :
        if log2fc_positive :
            df_result = df_result[ df_result[ label_log2fc ] > 0 ]
        else :
            df_result = df_result[ df_result[ label_log2fc ] < 0 ]
    return df_result


# In[ ]:


def GET_dataframe_Gene_Symbol_from_Gene_ID( df ) :
    ''' if index and columns are in same other, convert Gene_ID labels to Symbol labels and return the DataFrame '''
    if ( df.index.values == df.columns.values ).all( ) :
        Gene_Symbols = List_Gene_ID__2__List_Gene_Symbol( df.index.values ) 
        return pd.DataFrame( df.values, index = Gene_Symbols, columns = Gene_Symbols )
    else :
        return - 1


# ##### Functions for Pandas DataFrame containing Linear Data 

# In[ ]:


def DF_Tabular_2_Linear_with_filter( df, thres_lower = None, thres_higher = None ) :
    ''' Convert DataFrame containing Tabular data to linear data with filters for selecting data_values (upper and lower thresholds are available)  '''
    data = df.values
    arr_filter = np.ones_like( data ).astype( bool )
    if thres_lower is None and thres_higher is None :
        pass
    elif thres_lower is not None and thres_higher is not None :
        arr_filter = ( thres_higher < data ) | ( data < thres_lower )
    elif thres_lower is not None : 
        arr_filter = arr_filter & ( data < thres_lower )
    elif thres_higher is not None :
        arr_filter = arr_filter & ( thres_higher < data )
    s_indices, s_columns = pd.Series( DICTIONARY_Build_from_arr( df.index.values ) ), pd.Series( DICTIONARY_Build_from_arr( df.columns.values ) )
    indices, columns = np.where( arr_filter )
    data_values = data[ arr_filter ] 
    return pd.DataFrame( dict( Index = s_indices.loc[ indices ].values, Column = s_columns.loc[ columns ].values, Data_Value = data_values ) )#.sort_values( 'Correl_Coeffi', ascending = False )


# ### Functions for Pandas MultiIndex dataframe 

# In[ ]:


def PANDAS_MULTIINDEX_get_indices_from_multiIndex_of_df( df, index_among_multiIndex = 0 ) :
    ''' from a given multiIndexed dataframe, retrive and return i-th index of multiIndex as an numpy array (dtype = object)  '''
    arr_multiIndex = df.index.values
    if type( arr_multiIndex[ 0 ] ) is not tuple :
        return -1
    return np.array( list( multiIndex[ index_among_multiIndex ] for multiIndex in arr_multiIndex ), dtype = object )
MULTIINDEX_At = PANDAS_MULTIINDEX_get_indices_from_multiIndex_of_df


# In[ ]:


def PANDAS_MULTIINDEX_get_multiIndex_as_np_array( df ) :
    ''' from a given multiIndexed dataframe, retrive multiIndex tuple and return it as an numpy array (dtype = object)  '''
    arr_multiIndex = df.index.values
    if type( arr_multiIndex[ 0 ] ) is not tuple :
        return -1
    return np.array( list( list( index for index in multiIndex ) for multiIndex in arr_multiIndex ), dtype = object )


# ### Functions for Scipy.stats

# In[ ]:


def PEARSONR( arr_x, arr_y ) : # 2020-07-30 16:17:30 
    ''' remove np.nan values from two given arrays with the same length, and return the result of stats.pearsonr and number of pairs non-nan values. if number of non-nan values is smaller than 2, return np.nan, np.nan, 'n_non_nan_pairs' instead. '''
    mask_valid = ~ ( np.isnan( arr_x ) | np.isnan( arr_y ) )
    arr_x_valid, arr_y_valid = arr_x[ mask_valid ], arr_y[ mask_valid ] 
    n_non_nan_pairs = len( arr_x_valid )
    if n_non_nan_pairs > 1 : return ( * stats.pearsonr( arr_x_valid, arr_y_valid ), n_non_nan_pairs )
    else : return np.nan, np.nan, n_non_nan_pairs


# ### Normalization Functions

# ##### Normalization within entires

# In[ ]:


def ROW_NORMALIZE_log_mean( df ) :
    ''' Normalize each entry (row) by dividing log mean of data values of the entry '''
    df = deepcopy( df )
    df = df.replace( 0, np.nan ).dropna( )
    df.iloc[ :, : ] = ( df.values.T / np.power( 2, np.log2( df.values ).mean( axis = 1 ) ) ).T * 100 # normalize data by log-mean of each entry
    return df.dropna( )


# ##### Normalization between samples

# In[ ]:


def NORMALIZE_Relative_Log_Expression( df, sample_list_for_pseudoreference = None, sample_list = None, Gene_Set = None, allow_nan_values = False, use_median_ratio = True, scale_all_genes = True, dropna = False, return_normalized_ratios = False, normalization_description = '' ) :
    ''' Perform Relative Log Expression (RLE) Normalization on DataFrame ( col = sample, index = genes ) for sample_list (default : all columns ) and Gene_Set (default : all indices)
        if 'use_median_ratio' is False, use mean_ratio instead.
        return df and Series with median ratios to pseudoreference (log-averaged for 'sample_list_for_pseudoreference', default : all columns ) for the samples that has been set for normalization. 
        When 'return_normalized_ratios' is True, return normalized relative_ratios to pseudoreference instead of normalized data.
        When 'allow_nan_values' is set to True, use masked array to retrive scaling factors for each sample (if data contains a lot of NaN or zeros, it is better to set this argument to True) '''
    Gene_Set = LIST_intersection_with_set( Gene_Set, df.index.values ) if Gene_Set is not None else df.index.values # by default, Gene_Set is all genes in the df # if Gene_Set has been given, only consider genes that exist in df
    sample_list = LIST_intersection_with_set( sample_list, df.columns.values ) if sample_list is not None else df.columns.values # set default sample_list and sample_list_for_pseudoreference
    sample_list_for_pseudoreference = LIST_intersection_with_set( sample_list_for_pseudoreference, df.columns.values ) if sample_list_for_pseudoreference is not None else sample_list
    df = deepcopy( df ) # make a copy dataframe to be safe
    df_all_input_genes = df[ sample_list ].dropna( ) if dropna or return_normalized_ratios else df[ sample_list ] # drop NaN values for sample_list samples if 'dropna' is True or return_normalized_ratios
    df = df.loc[ Gene_Set, sample_list ].replace( 0, np.nan )
    df = df if allow_nan_values else df.dropna( ) # drop NaN and zero values for sample_list samples
    data = df.values.astype( float ) # row = gene, columns = samples
    data_pseudoref = df[ sample_list_for_pseudoreference ].values.astype( float )
    if allow_nan_values :
        data, data_pseudoref = np.ma.masked_array( data, np.isnan( data ) ), np.ma.masked_array( data_pseudoref, np.isnan( data_pseudoref ) )
    pseudoreference_values = np.power( 2, np.mean( np.log2( data_pseudoref ), axis = 1 ) ) # (pseudoreference) retrieve log-average values for each genes (by applying power of 1 / n_samples to a product of values of all samples)
    relative_ratios = ( data.T / pseudoreference_values ).T
    if use_median_ratio :
        scaling_ratios = np.ma.median( relative_ratios, axis = 0 ) if allow_nan_values else np.median( relative_ratios, axis = 0 ) # retrive median ratios of genes in Gene_Set to pseudoreference for each sample
    else :
        scaling_ratios = np.mean( relative_ratios, axis = 0 ) # retrive mean ratios of genes in Gene_Set to pseudoreference for each sample
    series_scaling_ratios = pd.Series( scaling_ratios.data if allow_nan_values else scaling_ratios, index = df.columns.values ) # build pandas series with median values with index as Sample_IDs 
    series_scaling_ratios = series_scaling_ratios.rename( 'Median '* use_median_ratio + 'scaling_ratios_used_for_normalization' + normalization_description ) # name a scale factor with a description given by 'normalization_description'
    if scale_all_genes : # if 'scale_all_genes' is set to True, return a data of all genes that are scaled with a median_ratio calculated from a given subset of genes 
        data, index, columns = df_all_input_genes.values, df_all_input_genes.index.values, df_all_input_genes.columns.values
        if return_normalized_ratios : # 'return_normalized_ratios' and 'scale_all_genes' are set to True, compute relative_ratios for all genes and scale them with median ratio
            data_pseudoref = df_all_input_genes[ sample_list_for_pseudoreference ].values.astype( float )
            if allow_nan_values :
                data_pseudoref = np.ma.masked_array( data_pseudoref, np.isnan( data_pseudoref ) )
            pseudoreference_values = np.power( 2, np.mean( np.log2( data_pseudoref ), axis = 1 ) )
            relative_ratios = ( data.T / pseudoreference_values ).T
    else :
        index, columns = df.index.values, df.columns.values
    data_norm = relative_ratios / scaling_ratios if return_normalized_ratios else data / scaling_ratios # divide relative ratios or data of each sample by its median ratio (scaling values)
    df_norm = pd.DataFrame( data_norm.data if allow_nan_values else data_norm, index = index, columns = columns )
    df_norm.index.name = df.index.name
    df_norm.columns.name = df.columns.name
    return df_norm, series_scaling_ratios # return median-ratio normalized data


# In[ ]:


def STANDARD_SCORE_Z_score( df, list_samples = None ):
    ''' Calculating z scores for list_samples using masked array and return dataframe with zsocres. Default list_samples = dict_type_sample[ 'all' ] '''
    if list_samples is None :
        list_samples = df.columns.values
    df_zscore = deepcopy ( df )
    data_values = df_zscore.loc[ :, list_samples ].values
    data_values = np.ma.masked_array( data_values, np.isnan( data_values ) )
    data_values_zscore = stats.zscore( data_values, axis = 1 )
    df_zscore.loc[ :, list_samples ] = data_values_zscore
    return df_zscore
Z_score = STANDARD_SCORE_Z_score


# ### Function for interacting with Proteome DataFrame

# In[ ]:


def GET_DATA_of_gene( Gene, df = None, sample_list = None ) :
    ''' retrive data of a given gene in a given df and sample_list. default df and sample_list are df_proteome and N_samples, respectively '''
    Gene_ID = Gene_2_Gene_ID( Gene ) # retrive Gene_ID of a given gene 
    if df is None : # set default df and sample_list 
        df = df_proteome
    if sample_list is None :
        sample_list = N_samples
    return np.array( df.loc[ Gene_ID, sample_list ].values, dtype = float ) # retrive data of a given gene in a given df and sample_list.


# ### Function for Accessing Gene_Sets of MSigDB Gene_Set Data

# In[ ]:


def GENE_SET_GET_Gene_Set_by_name( Gene_Set_name, return_Gene_Name_Symbol = False, gene_set_background = None ) :
    ''' return a Gene_Set (set of Gene_IDs) of a given Gene_Set_name in All Gene_Sets. if 'return_Gene_Name_Symbol' is True, return sorted list of Gene_Name_Symbol of a gene_set. If a 
    background gene_set has been given through 'gene_set_background' argument, return genes that exist in the background gene_set  '''
    for Gene_Sets_Name in dict_GeneSets_ALL.keys( ) :
        if Gene_Set_name in dict_GeneSets_ALL[ Gene_Sets_Name ] :
            Gene_Set = dict_GeneSets_ALL[ Gene_Sets_Name ][ Gene_Set_name ]
            if gene_set_background is not None : # If a background gene_set has been given through 'gene_set_background' argument, return genes that exist in the background gene_set
                Gene_Set = Gene_Set.intersection( gene_set_background )
            if return_Gene_Name_Symbol : # if 'return_Gene_Name_Symbol' is True, return sorted list of Gene_Name_Symbol of a gene_set.
                return List_Gene_ID__2__List_Gene_Symbol( Gene_Set, add_gene_name = True, return_sorted = True )
            else :
                return Gene_Set


# In[ ]:


"""def GENE_SET_GET_MSigDB_Gene_Set_by_name( Gene_Set_name, return_Gene_Name_Symbol = True, gene_set_background = None ) :
    ''' return a Gene_Set (set of Gene_IDs) of a given Gene_Set_name in MSigDB Gene_Sets. if 'return_Gene_Name_Symbol' is True, return sorted list of Gene_Name_Symbol of a gene_set. If a 
    background gene_set has been given through 'gene_set_background' argument, return genes that exist in the background gene_set  '''
    for Gene_Sets_Name in dict_MSigDB_name__set_name__gene_set.keys( ) :
        if Gene_Set_name in dict_MSigDB_name__set_name__gene_set[ Gene_Sets_Name ] :
            Gene_Set = dict_MSigDB_name__set_name__gene_set[ Gene_Sets_Name ][ Gene_Set_name ]
            if gene_set_background is not None : # If a background gene_set has been given through 'gene_set_background' argument, return genes that exist in the background gene_set
                Gene_Set = Gene_Set.intersection( gene_set_background )
            if return_Gene_Name_Symbol : # if 'return_Gene_Name_Symbol' is True, return sorted list of Gene_Name_Symbol of a gene_set.
                return List_Gene_ID__2__List_Gene_Symbol( Gene_Set, add_gene_name = True, return_sorted = True )
            else :
                return Gene_Set"""


# In[ ]:


def GENE_SET_Search_MSigDB_Gene_Set_Name( query ) :
    '''  Search Gene_Sets maching a given query in MSigDB GeneSets, and print the result, and return a list of GeneSet_names  '''
    Matched_Gene_Sets = list( ) # an empty list of matched gene_sets
    for GeneSets_name in dict_MSigDB_name__set_name__gene_set :
        Matched_Gene_Sets.extend( Search_list_of_strings( list( dict_MSigDB_name__set_name__gene_set[ GeneSets_name ].keys( ) ), query = query.upper( ) ) ) # retriving Gene_Sets matching the query # Since its MSiGDB genesets only, change query to upper characters 
    return Matched_Gene_Sets


# In[ ]:


def GENE_SET_Search_Gene_Set_Name_in_Gene_Sets( dict_Gene_Sets, query, return_GeneSet_names = True ) :
    '''  Search Gene_Sets maching a given query in a given dict_Gene_Sets, and print the result  '''
    Matched_Gene_Sets = Search_list_of_strings( list( dict_Gene_Sets.keys( ) ), query = query ) # retriving Gene_Sets matching the query
    print( 'number_of_genes', '\t', 'Gene_Set_Name' ) # print the number of genes for each matched gene set
    for Gene_Set in Matched_Gene_Sets :
        print( len( dict_Gene_Sets[ Gene_Set ] ), '\t', Gene_Set )
    return Matched_Gene_Sets


# In[ ]:


def GENE_SET_Search_Gene_in_Gene_Sets( dict_Gene_Sets, Gene, return_subset_dict_Gene_Sets = False ) :
    '''  Search for Gene_Sets having a query Gene in a given dict_Gene_Sets, and return the list of Gene_Set_Names as a list  '''
    Gene_ID = Gene_2_Gene_ID( Gene ) # retrive Gene_ID
    if Gene_ID == -1 : # if Gene is invalid, return an error value
        return -1
    list_Gene_Set_Names = list( Gene_Set_Name for Gene_Set_Name, Gene_Set in dict_Gene_Sets.items( ) if Gene_ID in Gene_Set )
    if return_subset_dict_Gene_Sets : # if 'return_subset_dict_Gene_Sets' is set to True, by using a pandas series indexing, return a subset of dict_Gene_Sets
        return pd.Series( dict_Gene_Sets ).loc[ list_Gene_Set_Names ].to_dict( ) 
    else :
        return list_Gene_Set_Names # return 


# In[ ]:


def GENE_SET_Search_Gene_in_Gene_Sets( dict_Gene_Sets, Gene ) :
    '''  Search for Gene_Sets having a query Gene in a given dict_Gene_Sets, and return the list of Gene_Set_Names as a list  '''
    Gene_ID = Gene_2_Gene_ID( Gene ) # retrive Gene_ID
    if Gene_ID == -1 : # if Gene is invalid, return an error value
        return -1
    return list( Gene_Set_Name for Gene_Set_Name, Gene_Set in dict_Gene_Sets.items( ) if Gene_ID in Gene_Set ) # return 


# In[ ]:


def GENE_SET_Get_a_union_of_GeneSets_from_list_of_MSigDB_GeneSet_names( list_MSigDB_GeneSet_names ) :
    ''' return a union (all unique elements) of MSigDB Gene_Sets from a given list of MSigDB_GeneSet_names  '''
    set_union = set( )
    for MSigDB_GeneSet_name in list_MSigDB_GeneSet_names :
        set_union.update( GENE_SET_GET_Gene_Set_by_name( MSigDB_GeneSet_name, return_Gene_Name_Symbol = False ) )
    return set_union


# ### Functions for GSEA analysis with KEGG pathway Gene_Sets and other Gene_Sets

# In[ ]:


def GSEA_pathway( pathway, Gene_Set, background_Gene_Set ) :
    '''
    Perform GSEA with KEGG pathway by using Fisher's exact test.
    '''
    if pathway not in dict_pathway :
        print( 'GSEA_pathway: invalid pathway name' )
        return -1
    num_gene_set = len( Gene_Set )
    num_background_genes = len( background_Gene_Set )
    pathway_gene_set = set( dict_pathway[ pathway ][ 'genes' ] )
    num_overlaped_genes = len( pathway_gene_set.intersection( Gene_Set ) )
    num_overlaped_background_genes = len( pathway_gene_set.intersection( background_Gene_Set ) )
    odd_ratio, p_value = fisher_exact( [ [ num_overlaped_background_genes, num_background_genes ], [ num_overlaped_genes, num_gene_set ] ] )
    return dict( num_overlaped_genes = num_overlaped_genes, num_overlaped_background_genes = num_overlaped_background_genes, odd_ratio = odd_ratio, p_value = p_value )


# In[ ]:


def KEGG_Pathways_GSEA( Gene_IDs, backgroud_genes ) :
    '''
    Perform GSEA on all KEGG pathway gene sets
    '''
    backgroud_gene_set = set( backgroud_genes )
    Gene_Set = set( Gene_IDs )
    dict_pathway_GSEA = { } # perform GSEA for every KEGG_pathway gene set and save result to a dictionary
    for pathway in dict_pathway : 
        if len( dict_pathway[ pathway ][ 'genes' ] ) != 0 :
            dict_pathway_GSEA[ pathway ] = GSEA_pathway( pathway, Gene_Set, backgroud_gene_set )
    df = pd.DataFrame( dict_pathway_GSEA ).transpose( )
    df = df.sort_values( [ 'p_value' ], axis = 'index' ) # sort result according to p_value
    # divide and a result dataframe into two by odd_ratio and return two dataframes inside a dictionary 
    return dict( over_represented = df[ df.odd_ratio < 1 ], under_represented = df[ df.odd_ratio > 1 ], result = df ) 


# In[ ]:


def GSEA_Gene_Set( Gene_Set, Gene_IDs, background_Gene_IDs ) :
    '''    Perform GSEA with a given Gene_Set by using Fisher's exact test.    '''
    num_genes = len( Gene_IDs )
    num_background_genes = len( background_Gene_IDs )
    Gene_Set = set( Gene_Set )
    num_overlaped_genes = len( Gene_Set.intersection( Gene_IDs ) )
    num_overlaped_background_genes = len( Gene_Set.intersection( background_Gene_IDs ) )
    odd_ratio, p_value = fisher_exact( [ [ num_overlaped_genes, num_genes ], [ num_overlaped_background_genes, num_background_genes ] ] )
    return dict( num_overlaped_genes = num_overlaped_genes, num_overlaped_background_genes = num_overlaped_background_genes, odd_ratio = odd_ratio, p_value = p_value )


# In[ ]:


def ANALYSIS_GSEA__dict_Gene_Sets( Gene_IDs, background_Gene_IDs = None, dict_Gene_Sets = None, thres_adj_p_val = 0.05 ) :
    '''    Perform GSEA using gene sets in a given dict_Gene_Sets and return results as two DataFrames
    Calculate adjusted p_value using fdr_bh method
    filter out entries that have higher adjusted_p_values than a given threshold adj_p_val    '''
    background_Gene_IDs = set( background_Gene_IDs )
    Gene_IDs = set( Gene_IDs )
    num_genes = len( Gene_IDs )
    num_background_genes = len( background_Gene_IDs )
    dict_GSEA = { } # perform GSEA for every gene set in a given dict_Gene_Sets and save result to a dictionary
    for Gene_Set_Name in dict_Gene_Sets : 
        Gene_Set = set( dict_Gene_Sets[ Gene_Set_Name ] )
        if len( Gene_Set ) == 0 :
            print( Gene_Set_Name, '; ERROR : length of gene set is zero' )
            return -1
        num_overlaped_genes = len( Gene_Set.intersection( Gene_IDs ) ) # perfoem GSEA of a current Gene_Set
        num_overlaped_background_genes = len( Gene_Set.intersection( background_Gene_IDs ) )
        odd_ratio, p_value = fisher_exact( [ [ num_overlaped_genes, num_genes ], [ num_overlaped_background_genes, num_background_genes ] ] )
        dict_GSEA[ Gene_Set_Name ] = dict( num_overlaped_genes = num_overlaped_genes, num_overlaped_background_genes = num_overlaped_background_genes, odd_ratio = odd_ratio, p_value = p_value )
    df = pd.DataFrame( dict_GSEA ).transpose( )
    df[ 'adj_p_val' ] = multipletests( df[ 'p_value' ].values , 0.05, 'fdr_bh', False, False )[ 1 ] # perform fdr correction of p_values
    df = df[ df.adj_p_val.values < thres_adj_p_val ] # filter out entries that have higher p_values than a given threshold p_value
    df = df.sort_values( [ 'p_value' ], axis = 'index' ) # sort result according to p_value
    return df # dict( over_represented = df[ df.odd_ratio > 1 ], under_represented = df[ df.odd_ratio < 1 ] ) # divide and a result dataframe into two by odd_ratio and return two dataframes inside a dictionary 


# In[ ]:


def ANALYSIS_GSEA__MSigDB_Gene_Sets( Gene_IDs, background_Gene_IDs = None, thres_adj_p_val = 0.05, Gene_Sets_Names = None ) :
    '''    Perform GSEA using all available MSigDB gene sets return results as DataFrames in a dictionary
    The list of Gene_Sets_Names can be given though an argument Gene_Sets_Names
    filter out entries that have higher adjusted_p_values than a given threshold adj_p_val    '''
    if background_Gene_IDs is None :
        background_Gene_IDs = PD_Select( df_ID_Symbol, Locus_Group = 'protein-coding gene' ).index.values # set protein-coding genes as a background geneset
    if Gene_Sets_Names is None : # if a list of Gene_Sets_Names was not given, set default Gene_Sets_Names
        Gene_Sets_Names = list( dict_MSigDB_name__set_name__gene_set.keys( ) )
    background_Gene_IDs = set( background_Gene_IDs )
    Gene_IDs = set( Gene_IDs )
    l_df = list( )
    for Gene_Sets_Name in Gene_Sets_Names : # for each Gene_Sets data in MSigDB datasets
        print( Gene_Sets_Name ) # report the name of Gene_Sets data that will be used for GSEA 
        df_result = ANALYSIS_GSEA__dict_Gene_Sets( Gene_IDs = Gene_IDs, background_Gene_IDs = background_Gene_IDs, dict_Gene_Sets = dict_MSigDB_name__set_name__gene_set[ Gene_Sets_Name ], thres_adj_p_val = thres_adj_p_val )
        if isinstance( df_result, pd.DataFrame ) :
            df_result[ 'Gene_Sets_Name' ] = Gene_Sets_Name
            l_df.append( df_result )
    df = pd.concat( l_df ) # combine results from multiple gene_sets 
    df[ 'adj_p_val' ] = multipletests( df[ 'p_value' ].values , 0.05, 'fdr_bh', False, False )[ 1 ] # perform multiple correction again using p_values 
    return df # return a dictionary containing result DataFrames


# ### Statistic Utility functions

# In[ ]:


def MULTIPLE_TEST_p_value( arr_p_values ) :
    ''' adjust p_values for only valid p_values in given array of arr_p_values. NaN values are considered invalid p_values '''
    arr_adj_p_values = np.ones_like( arr_p_values ) 
    mask_valid_p_value = ~ np.isnan( arr_p_values ) 
    arr_adj_p_values[ mask_valid_p_value ] = multipletests( arr_p_values[ mask_valid_p_value ], 0.05, 'fdr_bh', False, False )[ 1 ]
    return arr_adj_p_values


# ### Calculate Changes in Tumor compare to Normal 

# In[ ]:


def Change__normal_vs_tumor__Gene_Set( df, Gene_Set = None, print_error_msg = False ) :
    '''
    Calculate of Log2FC (ratio of tumor to normal) of averages of protein amounts of a given set of genes, 
    Categoties are 'normal' tissue and 'tumor' tissue
    '''
    # process input gene set 
    if Gene_Set is not None :
        Gene_List = list( Gene_Set )
        num_genes =  len( Gene_List )
    else :
        print( 'invalid inputs' )
        return -1
    if num_genes == 0 :
        print( 'gene_set is empty (invalid gene_set)')
        return -1
    # select genes existing in proteome data
    existing_gene_set = set( df.index.values ).intersection( Gene_List )
    # calculate number of valid genes and total number of genes of the pathway or set and print the values 
    valid_gene_set = set( df[ dict_type_sample[ 'normal_or_tumor' ] ].loc[ existing_gene_set ].dropna( ).index ) # only process genes that do not contain NaN values 
    num_valid_genes = len( valid_gene_set )
    if num_valid_genes == 0 : # if there is no valid genes, end the method
        if print_error_msg :
            print( 'number of valid genes = 0' )
        return -1 
    coverage = round( num_valid_genes / num_genes * 100, 1 ) # calculate pathway or gene set coverage
    # data values 
    # compute average of values to assess Gene_Set activity
    normal_data = np.average( df.loc[ valid_gene_set, dict_type_sample[ 'normal' ] ].values, axis = 0 )
    tumor_data = np.average( df.loc[ valid_gene_set, dict_type_sample[ 'tumor' ] ].values, axis = 0 )
    # perform T test
    p_value = stats.ttest_ind( normal_data, tumor_data )[ 1 ]
    Log2FC_ratio_of_tumor_to_normal = np.log2( np.average( tumor_data ) / np.average( normal_data ) )
    return dict( num_valid_genes = num_valid_genes, num_genes = num_genes, p_value = p_value, 
                 Log2FC_ratio_of_tumor_to_normal = Log2FC_ratio_of_tumor_to_normal, coverage = coverage )


# In[ ]:


def Calculate_change__dict_Gene_Sets( df, dict_Gene_Sets, print_error_msg = False ) :
    '''    Calculate activities of Gene_Sets     '''
    dict_Gene_Sets__normal_vs_tumor = { }
    for Gene_Set_name in dict_Gene_Sets :
        Gene_Set = dict_Gene_Sets[ Gene_Set_name ]
        if len( Gene_Set ) != 0 : # if current gene_set is valid
            result = Change__normal_vs_tumor__Gene_Set( df, Gene_Set = Gene_Set, print_error_msg = print_error_msg )
            if result != -1 : # check whether there has been an error
                dict_Gene_Sets__normal_vs_tumor[ Gene_Set_name ] = result
    df = pd.DataFrame( dict_Gene_Sets__normal_vs_tumor ).transpose()
    # FDR correction # perform one of the most strict FDR corrections (since the number of multiple testing is not large)
    p_values = np.array( df[ 'p_value' ].values, dtype = float )
    df[ 'adjuested_p_value' ] = multipletests( p_values , 0.05, 'bonferroni', False, False)[ 1 ] 
    # sort and convert it to excel file
    df = df.sort_values( [ 'Log2FC_ratio_of_tumor_to_normal' ], axis = 'index' )
    return df


# In[ ]:


def Calculate_min_by_max( df, sample_list, Gene_Set = None, n_std_for_outliers = 3 ) :
    ''' Calculate ratio of min to max values of a given data ('sample_list' in 'df') for genes in 'Gene_Set' while excluding outliers outside 'n_std_for_outliers'. Return result as a series  '''
    if Gene_Set is not None :
        df = PANDAS_Subset( df, Gene_Set )
    df = df[ sample_list ] 
    data = df.values
    data = np.ma.masked_array( data, OUTLIERS_GET_mask_for_outliers( data, n_std_for_outliers = n_std_for_outliers ) )
    return pd.Series( data.min( axis = 1 ) / data.max( axis = 1 ), index = df.index.values )


# In[ ]:


def CALCULATE_masked_arr__average__std_by_mean( arr ) :
    ''' For each row, mask NaN values, calculate average and std / average using the masked array. return arr_masked, arr_average, arr_std_by_mean '''
    if arr.ndim == 1 :
        arr = np.ma.masked_array( arr, np.isnan( arr ) ) # mask their NaN values for downstream calculations
        return arr, arr, np.zeros_like( arr )
    else :
        arr = np.ma.masked_array( arr, np.isnan( arr ) ) # mask their NaN values for downstream calculations
        arr_average = np.average( arr, axis = 1 )
        arr_std_by_mean = np.std( arr, axis = 1 ) / arr_average # calculate std divided by mean
        return arr, arr_average, arr_std_by_mean


# In[ ]:


def Calculate_Log2FC_p_value__A_vs_B( df, Gene_IDs = None, sample_list_condition_A = None, sample_list_condition_B = None, p_value_ttest_rel = False, label_condition_A = 'condition_A', label_condition_B = 'condition_B', method_average = 'log_mean', add_raw_data = False, replace_zero_with = None, drop_avg = False, drop_std = True ) :
    ''' Calculate Log2FC, t-test p_values, and fdr-bh adjusted p_values for condition 'A' (usually Normal) and condition 'B' (usually Tumor)
        samples for a given Gene_IDs and return the result as a dataframe. Log2FC is log2( condition B average / condition A average )
        use masked array to avoid propagation of NaN values to Log2FC and p_values. label columns according to given labels 'label_condition_A' and 'label_condition_B'.
        If 'replace_zero_with' is not None, replace zero values in df (Gene_IDs subset of df) with a value given by 'replace_zero_with' and drop NaN values. Add original data if 'add_raw_data' is True 
        'method_average' = 'log_mean' or 'mean' '''
    Gene_IDs = df.index.values if Gene_IDs is None else LIST_intersection_with_set( Gene_IDs, df.index.values )  # if Gene_IDs was not given, set indices of df as Gene_IDs
    if sample_list_condition_A is None :
        sample_list_condition_A = dict_type_sample[ 'normal' ]
    if sample_list_condition_B is None :
        sample_list_condition_B = dict_type_sample[ 'tumor' ]
    set_all_samples = set( df.columns.values )
    sample_list_condition_A, sample_list_condition_B = LIST_intersection_with_set( sample_list_condition_A, set_all_samples ), LIST_intersection_with_set( sample_list_condition_B, set_all_samples ) # retrive only valid samples
    print( "Number of valid sample A: {}, sample B: {}".format( len( sample_list_condition_A ), len( sample_list_condition_B ) ) ) # print number of valid samples
    if len( sample_list_condition_A ) != len( sample_list_condition_B ) : # if valid numbers of samples of condition A and B are not the same, calculate ttest_independent using ttest_ind function
        p_value_ttest_rel = False
    if len( sample_list_condition_A ) == 0 or len( sample_list_condition_B ) == 0 : # if number of samples of one condition is zero, return -1
        return -1
    all_samples = list( sample_list_condition_A ) + list( sample_list_condition_B ) # create a sample list with all samples
    df = df.loc[ Gene_IDs ] # retrive data of a given Gene_IDs
    if replace_zero_with is not None :
        df = df.replace( 0, replace_zero_with ).dropna( )
    Gene_IDs = df.index.values # retrieve Gene_IDs from the given dataframe
    condition_A_data, condition_B_data, data_all = df[ sample_list_condition_A ].values,  df[ sample_list_condition_B ].values, df[ all_samples ].values
    condition_A_data, condition_B_data = np.ma.masked_array( condition_A_data, np.isnan( condition_A_data ) ), np.ma.masked_array( condition_B_data, np.isnan( condition_B_data ) ) # mask their NaN values for downstream calculations
    data_all = np.ma.masked_array( data_all, np.isnan( data_all ) )
    if method_average == 'mean' :
        condition_A_average, condition_B_average, average_all = np.average( condition_A_data, axis = 1 ), np.average( condition_B_data, axis = 1 ), np.average( data_all, axis = 1 )
    else :
        condition_A_average, condition_B_average, average_all = np.power( 2, np.log2( condition_A_data ).mean( axis = 1 ) ), np.power( 2, np.log2( condition_B_data ).mean( axis = 1 ) ), np.power( 2, np.log2( data_all ).mean( axis = 1 ) )
    condition_A_std_by_mean, condition_B_std_by_mean = np.std( condition_A_data, axis = 1 ) / condition_A_average, np.std( condition_B_data, axis = 1 ) / condition_B_average # calculate std divided by mean
    std_by_mean_all = np.std( data_all, axis = 1 ) / average_all
    if len( sample_list_condition_A ) == 1 or len( sample_list_condition_B ) == 1 : # if number of samples of one condition is zero, p_values and adj_p_values are np.nan
        p_values, adj_p_values = np.full( len( df ), np.nan ), np.full( len( df ), np.nan )
    else :
        p_values = mstats.ttest_rel( condition_A_data, condition_B_data, axis = 1 ).pvalue if p_value_ttest_rel else mstats.ttest_ind( condition_A_data, condition_B_data, axis = 1 ).pvalue # perform the t-test (row by row) # use mstats module which supports masked arrays
        adj_p_values = MULTIPLE_TEST_p_value( p_values ) # calculate adjusted p_values
    if label_condition_A == 'cl name' or label_condition_B == 'cl name' :
        label_condition_A, label_condition_B = '_'.join( list( cl.split( '_' )[ 0 ] for cl in sample_list_condition_A ) ), '_'.join( list( cl.split( '_' )[ 0 ] for cl in sample_list_condition_B ) )
    dict_rename = dict( condition_A_average = '{A}_average'.format( A = label_condition_A ), condition_B_average = '{B}_average'.format( B = label_condition_B ), condition_A_std_by_mean = '{A}_std_by_mean'.format( A = label_condition_A ), condition_B_std_by_mean = '{B}_std_by_mean'.format( B = label_condition_B ) ) # define dictionary for rename columns according to given condition A and B labels
    Log2FC = condition_B_average / condition_A_average
    Log2FC = np.ma.masked_array( Log2FC, mask = np.isnan( Log2FC ) )
    Log2FC = np.ma.log2( Log2FC ) # calculate log2 fold change values
    Difference = condition_B_average - condition_A_average # calculate difference
    df_result = pd.DataFrame( dict( Log2_Fold_Change = Log2FC, Difference = Difference, p_value = p_values, adjusted_p_value = adj_p_values, average_all = average_all, condition_A_average = condition_A_average, condition_B_average = condition_B_average, condition_A_std_by_mean = condition_A_std_by_mean, condition_B_std_by_mean = condition_B_std_by_mean, std_by_mean_all = std_by_mean_all ), index = Gene_IDs )
    if drop_avg : # drop unneccesary columns
        df_result = df_result.drop( columns = [ 'average_all', 'condition_A_average', 'condition_B_average' ] )
    if drop_std :
        df_result = df_result.drop( columns = [ 'condition_A_std_by_mean', 'condition_B_std_by_mean', 'std_by_mean_all' ] )
    if len( sample_list_condition_A ) == 1 or len( sample_list_condition_B ) == 1 : # if number of samples of one condition is zero, drop p_value columns
        df_result = df_result.drop( columns = [ 'p_value', 'adjusted_p_value' ] )
    df_result = df_result.rename( columns = dict_rename ).sort_values( 'Log2_Fold_Change' ) # build a sorted result dataframe with labeled columns 
    if add_raw_data :
        df_result = df_result.join( df[ all_samples ] )
    if type( Gene_IDs[ 0 ] ) is tuple : # if index is a multi-index (tuple type when returned by df.index.values), convert index of df into MultiIndex
        print( 'Multi-Indexed df was received, return Multi-indexed df_result' ) 
        df_result.index = pd.MultiIndex.from_tuples( df_result.index )
        return df_result
    else :
        return df_result


# In[ ]:


def Calculate_Log2FC_p_value__dict_Gene_Sets__A_vs_B( df, dict_Gene_Sets, Gene_IDs = None, sample_list_condition_A = None, sample_list_condition_B = None ) :
    ''' Calculate Log2FC, t-test p_values, and fdr-bh adjusted p_values for condition 'A' (default: Normal) and condition 'B' (default: Tumor)
        samples for Gene_Sets in 'dict_Gene_Sets' with a given Gene_IDs and return the result as a dataframe. Log2FC is log2( condition B average / condition A average )
        dropna values to allow pair-wise comparison and calculation of Log2FC of genes in a Gene_Set    '''
    if Gene_IDs is None : # if Gene_IDs was not given, set indices of df as Gene_IDs
        Gene_IDs = df.index.values
    else : # if Gene_IDs has been given
        Gene_IDs = set( df.index.values ).intersection( Gene_IDs ) # consider Gene_IDs that exist in df
    if sample_list_condition_A is None :
        sample_list_condition_A = dict_type_sample[ 'normal' ]
    if sample_list_condition_B is None :
        sample_list_condition_B = dict_type_sample[ 'tumor' ]
    all_samples = list( sample_list_condition_A ) + list( sample_list_condition_B ) # create a sample list with all samples
    All_valid_Gene_IDs = set( df.loc[ Gene_IDs ].dropna( ).index.values ) # retrieve Gene_IDs without NaN from the given dataframe
    df = df.loc[ All_valid_Gene_IDs ] # retrive data of a given Gene_IDs
    dict_Gene_Set_result = dict( )
    for Gene_Set_name in dict_Gene_Sets : # for each Gene_Set_name
        Gene_Set = All_valid_Gene_IDs.intersection( dict_Gene_Sets[ Gene_Set_name ] ) # retrive valid Gene_IDs
        if len( Gene_Set ) == 0 : # if number of valid Genes is zero, skip the gene_set
            continue 
        df_Gene_Set = df.loc[ Gene_Set ] # retrive data of Gene_Set
        condition_A_data, condition_B_data, data_all = df_Gene_Set[ sample_list_condition_A ].values,  df_Gene_Set[ sample_list_condition_B ].values, df_Gene_Set[ all_samples ].values
        condition_A_average, condition_B_average, average_all = np.average( condition_A_data, axis = 1 ), np.average( condition_B_data, axis = 1 ), np.average( data_all, axis = 1 ) # calculate average of condition A samples and condition B samples and all samples for each gene
        condition_A_std_by_mean, condition_B_std_by_mean = ( np.std( condition_A_data, axis = 1 ) / condition_A_average ).mean( ), ( np.std( condition_B_data, axis = 1 ) / condition_B_average ).mean( ) # calculate std in condition A and condition B
        coondition_A_B_std_by_mean_avg = ( condition_A_std_by_mean + condition_B_std_by_mean ) / 2 # average std_by_mean in condition A and std in condition B, which can be used in an interactive plotting
        Log2FC = np.log2( condition_B_average / condition_A_average ) # calculate log2 fold change values
        p_value = stats.ttest_rel( condition_A_average, condition_B_average ).pvalue # perform the t-test (pairwise t-test, t-test of related samples)
        dict_Gene_Set_result[ Gene_Set_name ] = dict( Log2_Fold_Change = Log2FC.mean( ), Log2_Fold_Change_std = Log2FC.std( ), p_value = p_value, condition_A_average = condition_A_average.mean( ), condition_B_average = condition_B_average.mean( ), number_valid_genes = len( Gene_Set ), condition_A_std_by_mean = condition_A_std_by_mean, condition_B_std_by_mean = condition_B_std_by_mean, coondition_A_B_std_by_mean_avg = coondition_A_B_std_by_mean_avg )
    df_result = pd.DataFrame( dict_Gene_Set_result ).T # convert result to DataFrame
    df_result[ 'adjusted_p_value' ] = MULTIPLE_TEST_p_value( df_result.p_value.values ) # calculate adjusted p_values
    return df_result.sort_values( 'p_value' ) # sort results according to p_values and return result dataframe        


# In[ ]:


def Calculate_Log2FC_p_value___dict_Gene_Sets_Name__Gene_Sets___A_vs_B( df, dict_Gene_Sets_Name__Gene_Sets = None, Gene_IDs = None, filter_b4_cal__p_value = 0.05, filter_b4_cal__avg_all = 0, sample_list_condition_A = None, sample_list_condition_B = None, label_condition_A = 'condition_A', label_condition_B = 'condition_B', add_Gene_Set_2_df = True ) :
    ''' Calculate Log2FC, t-test p_values, and fdr-bh adjusted p_values for condition 'A' (default: Normal) and condition 'B' (default: Tumor)
        samples for each Gene_Sets in 'dict_Gene_Sets_Name__Gene_Sets' (default : 'dict_GeneSets_ALL' ) with a given Gene_IDs and return the result as a dataframe. 
        Log2FC is log2( condition B average / condition A average ). dropna values to allow pair-wise comparison and calculation of Log2FC of genes in a Gene_Set.
        label columns according to given labels 'label_condition_A' and 'label_condition_B'.'''
    df = deepcopy( df ) # copy dataframe to avoid changing original content
    if Gene_IDs is None : # if Gene_IDs was not given, set indices of df as Gene_IDs
        Gene_IDs = df.index.values
    else : # if Gene_IDs has been given 
        Gene_IDs = set( df.index.values ).intersection( Gene_IDs ) # consider Gene_IDs that exist in df
    if sample_list_condition_A is None :
        sample_list_condition_A = dict_type_sample[ 'normal' ]
    if sample_list_condition_B is None :
        sample_list_condition_B = dict_type_sample[ 'tumor' ]
    if dict_Gene_Sets_Name__Gene_Sets is None :
        dict_Gene_Sets_Name__Gene_Sets = dict_GeneSets_ALL
    set_all_samples = set( df.columns.values )
    sample_list_condition_A, sample_list_condition_B = list( set_all_samples.intersection( sample_list_condition_A ) ), list( set_all_samples.intersection( sample_list_condition_B ) ) # retrive only valid samples
    all_samples = list( sample_list_condition_A ) + list( sample_list_condition_B ) # create a sample list with all samples
    if ( len( sample_list_condition_A ) > 1 and len( sample_list_condition_B ) > 1 ) and ( filter_b4_cal__p_value is not None or filter_b4_cal__avg_all is not None ) : # if number of samples for each condition > 1, filter entries of df by p_value of Log2FC before Gene-Set based calculation, since the analysis assumes Log2FC of the genes to be significant
        df_result = Calculate_Log2FC_p_value__A_vs_B( df, sample_list_condition_A = sample_list_condition_A, sample_list_condition_B = sample_list_condition_B, drop_avg = False )
        df = df.loc[ RESULT_filter( df_result, thres_p_value = filter_b4_cal__p_value, thres_avg_all = filter_b4_cal__avg_all ).index.values ]
        print( "Number of Remaining Entries after filter :", len( df ) )
    All_valid_Gene_IDs = set( PD_Subset( df, Gene_IDs ).dropna( ).index.values ) # retrieve Gene_IDs without NaN from the given dataframe
    df = df.loc[ All_valid_Gene_IDs ] # retrive data of a given Gene_IDs
    ma_data = np.ma.masked_array( df.values, df.values == 0 ) # mask values with zero data values
    df.iloc[ :, : ] = ( ma_data.T / np.power( 2, np.log2( ma_data ).mean( axis = 1 ) ) ).T * 100 # normalize data by log-mean of each entry
    list_df_result = list( )
    for Gene_Sets_name, dict_Gene_Sets in dict_Gene_Sets_Name__Gene_Sets.items( ) :
        dict_Gene_Set_result = dict( )
        for Gene_Set_name in dict_Gene_Sets : # for each Gene_Set_name
            Gene_Set = All_valid_Gene_IDs.intersection( dict_Gene_Sets[ Gene_Set_name ] ) # retrive valid Gene_IDs
            if len( Gene_Set ) == 0 : # if number of valid Genes is zero, skip the gene_set
                continue 
            df_Gene_Set = PD_Subset( df, Gene_Set ) # retrive data of Gene_Set
            condition_A_data, condition_B_data, data_all = df_Gene_Set[ sample_list_condition_A ].values,  df_Gene_Set[ sample_list_condition_B ].values, df_Gene_Set[ all_samples ].values
            condition_A_average, condition_B_average, average_all = np.average( condition_A_data, axis = 1 ), np.average( condition_B_data, axis = 1 ), np.average( data_all, axis = 1 ) # calculate average of condition A samples and condition B samples and all samples for each gene
            condition_A_std_by_mean_all_genes, condition_B_std_by_mean_all_genes = ( np.std( condition_A_data, axis = 1 ) / condition_A_average ), ( np.std( condition_B_data, axis = 1 ) / condition_B_average )  # calculate std in condition A and condition B
            condition_A_std_by_mean, condition_B_std_by_mean = condition_A_std_by_mean_all_genes[ condition_A_average != 0 ].mean( ), condition_B_std_by_mean_all_genes[ condition_B_average != 0 ].mean( ) # calculate average of std/mean for all genes (except when average values across samples are zero)
            all_std_by_mean = ( np.std( data_all, axis = 1 ) / average_all ).mean( ) # average std_by_mean in condition A and std in condition B, which can be used in an interactive plotting
            log2FC = np.log2( condition_B_average / condition_A_average ) # calculate log2 fold change values
            log2FC = np.ma.masked_array( log2FC, np.isinf( log2FC ) ) # mask infinite values in log2FC (caused by 0 count values)
            p_value = stats.ttest_rel( condition_A_average, condition_B_average ).pvalue # perform the t-test (pairwise t-test, t-test of related samples)
            dict_Gene_Set_result[ Gene_Set_name ] = dict( Log2_Fold_Change = log2FC.mean( ), Log2_Fold_Change_std = log2FC.std( ), p_value = p_value, condition_A_average = condition_A_average.mean( ), condition_B_average = condition_B_average.mean( ), number_valid_genes = len( df_Gene_Set ), condition_A_std_by_mean = condition_A_std_by_mean, condition_B_std_by_mean = condition_B_std_by_mean, all_std_by_mean = all_std_by_mean )
            if add_Gene_Set_2_df : # if 'add_Gene_Set_2_df' is set True, add used Gene_Set to df for subsequent use 
                dict_Gene_Set_result[ Gene_Set_name ][  'Gene_Set' ] = Gene_Set
        df_result = pd.DataFrame( dict_Gene_Set_result ).T # convert result to DataFrame
        df_result[ 'Gene_Sets_Name' ] = np.full( len( df_result ), Gene_Sets_name )
        list_df_result.append( df_result )
    df_result = list_df_result[ 0 ] # if only one Gene_Sets is given, take the df_result out of the one-element list
    if len( list_df_result ) > 1 : # if more than one Gene_Sets were given, join the results in list_df_result into one DataFrame
        for a_df in list_df_result[ 1: ] :
            df_result = df_result.T.join( a_df.T ).T
    dict_rename = dict( condition_A_average = '{A}_average'.format( A = label_condition_A ), condition_B_average = '{B}_average'.format( B = label_condition_B ), condition_A_std_by_mean = '{A}_std_by_mean'.format( A = label_condition_A ), condition_B_std_by_mean = '{B}_std_by_mean'.format( B = label_condition_B ) ) # define dictionary for rename columns according to given condition A and B labels
    df_result[ 'adjusted_p_value' ] = MULTIPLE_TEST_p_value( df_result.p_value.values.astype( float ) ) # calculate adjusted p_values
    return df_result.sort_values( 'p_value' ).rename( columns = dict_rename ), df # sort results according to p_values and return result dataframe        


# In[ ]:


def GSEA_Calculate_Average___dict_Gene_Sets( series, dict_Gene_Sets, Gene_IDs = None, log10_transform = True, name_of_series = 'data' ) :
    ''' if name of a given series is not None, use the name to label columns, or, use name_of_series to label columns  '''
    if Gene_IDs is None : # if Gene_IDs was not given, set indices of df as Gene_IDs
        Gene_IDs = series.index.values
    else : # if Gene_IDs has been given
        Gene_IDs = set( series.index.values ).intersection( Gene_IDs ) # consider Gene_IDs that exist in df
    Gene_IDs = set( series.loc[ Gene_IDs ].dropna( ).index.values ) # retrieve Gene_IDs without NaN from the given dataframe
    series = series.loc[ Gene_IDs ] # retrive data of a given Gene_IDs
    dict_Gene_Set_result = dict( )
    for Gene_Set_name in dict_Gene_Sets : # for each Gene_Set_name
        Gene_Set = Gene_IDs.intersection( dict_Gene_Sets[ Gene_Set_name ] ) # retrive valid Gene_IDs
        if len( Gene_Set ) == 0 : # if number of valid Genes is zero, skip the gene_set
            continue 
        data = series.loc[ Gene_Set ]
        if log10_transform :
            data = np.log10( data )
        dict_Gene_Set_result[ Gene_Set_name ] = dict( Average = data.mean( ), Standard_dev = data.std( ), number_valid_genes = len( data ) )
    if series.name is not None :
        name_of_series = series.name
    Average_label, Std_label = 'Log10_' * log10_transform + 'Average_' + name_of_series, 'Standard_dev_of_' + 'Log10_' * log10_transform + name_of_series
    df_result = pd.DataFrame( dict_Gene_Set_result ).T.rename( columns = dict( Average = Average_label, Standard_dev = Std_label ) ) # convert result to DataFrame
    return df_result.sort_values( Average_label ) # sort results according to p_values and return result dataframe        


# In[ ]:


def GSEA_Calculate_Average___dict_Gene_Sets_Name__Gene_Sets( series, dict_Gene_Sets_Name__Gene_Sets = None, Gene_IDs = None, log10_transform = False, name_of_series = 'data' ) :
    ''' 'method_for_average' : 'average' or 'log_average' or a method   '''
    if dict_Gene_Sets_Name__Gene_Sets is None :
        dict_Gene_Sets_Name__Gene_Sets = dict_MSigDB_name__set_name__gene_set
    if Gene_IDs is None : # if Gene_IDs was not given, set indices of df as Gene_IDs
        Gene_IDs = series.index.values
    else : # if Gene_IDs has been given
        Gene_IDs = set( series.index.values ).intersection( Gene_IDs ) # consider Gene_IDs that exist in df
    Gene_IDs = set( series.loc[ Gene_IDs ].dropna( ).index.values ) # retrieve Gene_IDs without NaN from the given dataframe
    series = series.loc[ Gene_IDs ] # retrive data of a given Gene_IDs
    if series.name is not None : # use name of the series to label the columns of result dataframe 
        name_of_series = series.name
    Average_label, Std_label = 'Log10_' * log10_transform + 'Average_' + name_of_series, 'Standard_dev_of_' + 'Log10_' * log10_transform + name_of_series # set lables 
    list_df_result = list( )
    for Gene_Sets_name, dict_Gene_Sets in dict_Gene_Sets_Name__Gene_Sets.items( ) :
        dict_Gene_Set_result = dict( )
        for Gene_Set_name in dict_Gene_Sets : # for each Gene_Set_name
            Gene_Set = Gene_IDs.intersection( dict_Gene_Sets[ Gene_Set_name ] ) # retrive valid Gene_IDs
            if len( Gene_Set ) == 0 : # if number of valid Genes is zero, skip the gene_set
                continue 
            data = series.loc[ Gene_Set ]
            if log10_transform :
                data = np.log10( data )
            dict_Gene_Set_result[ Gene_Set_name ] = dict( Average = data.mean( ), Standard_dev = data.std( ), number_valid_genes = len( data ), Gene_Set = Gene_Set )
        df_result = pd.DataFrame( dict_Gene_Set_result ).T.rename( columns = dict( Average = Average_label, Standard_dev = Std_label ) ) # convert result to DataFrame
        df_result[ 'Gene_Sets_Name' ] = np.full( len( df_result ), Gene_Sets_name )
        list_df_result.append( df_result )
    df_result = list_df_result[ 0 ] # if only one Gene_Sets is given, take the df_result out of the one-element list
    if len( list_df_result ) > 1 : # if more than one Gene_Sets were given, join the results in list_df_result into one DataFrame
        for a_df in list_df_result[ 1 : ] :
            df_result = df_result.T.join( a_df.T ).T
    return df_result.sort_values( Average_label ) # sort results according to p_values and return result dataframe      


# ### Functions using sklearn modules

# In[ ]:


def CALCULATE_Factor_Analysis__df( df, sample_list, n_factors, Gene_Set = None ) :
    ''' using the given df and sample_list, drop NaN values and perform Factor Analysis and return DataFrame with the result '''
    if Gene_Set is None : # set default Gene_Set
        Gene_Set = df.index.values
    df = PANDAS_Subset( df, Gene_Set ) # subset df for Gene_Set
    df = df[ sample_list ].dropna( )
    column_labels = 'Factor_' + np.array( np.array( np.arange( 1, n_factors + 1 ), dtype = str ), dtype = object ) 
    FA = FactorAnalysis( n_components = n_factors )
    return pd.DataFrame( FA.fit_transform( df.values ) , index = df.index.values, columns = column_labels )


# ### Functions for Assisting Plotting Functions 

# In[ ]:


def GENE_Annotate_df_for_plot( df, gene_list_annotation ) :
    '''  gene_list_annotation : if a list of Gene is given through this argument and df is not a list (a single DataFrame), draw annotation of valid genes on the plot. If gene_list_annotation = 'all', 
    annotate all genes in a given df. If 'gene_list_annotation' is integer, annotate every 'gene_list_annotation' number of genes in df. if 'gene_list_annotation' is list containing three
    integers (the list should have at least two integers), annotate [ gene_list_annotation[ 0 ] : gene_list_annotation[ 1 ] : gene_list_annotation[ 2 ] (default : 1) ] genes in df. '''
    if gene_list_annotation is not None : # if 'gene_list_annotation' is given, retrive valid Gene_IDs and Gene_Symbols from the given list, and retrive positions of the valid genes in the DataFrame  
        Gene_IDs, list_positions = df.index.values, np.arange( 0, len( df.index.values ) ) # retrive all Gene_IDs and list of positions
        dict_ID_positions = dict( ( Gene_ID, position ) for Gene_ID, position in zip( Gene_IDs, list_positions ) ) # build dictionary ( key = Gene_ID, value = position ) 
        if gene_list_annotation == 'all' : # If gene_list_annotation = 'all', annotat all genes in a given df
            gene_list_annotation = Gene_IDs 
        elif type( gene_list_annotation ) is int :
            gene_list_annotation = Gene_IDs[ : : gene_list_annotation ]
        elif type( gene_list_annotation ) is list and type( gene_list_annotation[ 0 ] ) is int and len( gene_list_annotation ) < 4 : 
            if len( gene_list_annotation ) == 2 :
                gene_list_annotation.append( 1 ) # by default gene_list_annotation[ 2 ] = 1
            gene_list_annotation = Gene_IDs[ gene_list_annotation[ 0 ] : gene_list_annotation[ 1 ] : gene_list_annotation[ 2 ] ]
        if isinstance( Gene_IDs[ 0 ], ( int, float, np.float64, np.int64 ) ) : # if a given dataframe contains Entrez Gene_IDs
            Gene_IDs_annotation = LIST_intersection_with_set( List_Gene__2__List_Gene_ID( gene_list_annotation ), Gene_IDs ) # retrive valid Gene_IDs in a given list of genes for annotation
            Symbols_annotation = List_Gene_ID__2__List_Gene_Symbol( Gene_IDs_annotation ) # retrive valid Gene_Symbols from Gene_IDs
        else : # if data type of ID of entries is not integer or float (not Entrez Gene ID), annotate entries by using ID of entries directly
            Gene_IDs_annotation = LIST_intersection_with_set( gene_list_annotation, Gene_IDs )
            Symbols_annotation = Gene_IDs_annotation
        list_positions_annotation = list( dict_ID_positions[ Gene_ID ] for Gene_ID in Gene_IDs_annotation ) # retrive positions of the valid genes in the DataFrame  
        return Symbols_annotation, list_positions_annotation
    else :
        return -1, -1


# ### Functions for Clustering Correlation Matrix and its Visualization

# In[ ]:


def fancy_dendrogram( *args, truncate_mode = 'lastp', p = 12, leaf_rotation = 90., leaf_font_size = 12., show_contracted = True, annotate_above = 10, **kwargs ):
    max_d = kwargs.pop('max_d', None) # put defulat keyworded argument value
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = sch.dendrogram(*args, truncate_mode = truncate_mode, p = p, leaf_rotation = leaf_rotation, leaf_font_size = leaf_font_size, show_contracted = show_contracted, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


# In[ ]:


def PLOT_elbow__find_optimal_n_cluster_from_linkage_map( L, n_max_clusters = 10 ) :
    '''     Plot distance from root for last n_max_clusters number of clusters and 2nd derivative of it, 
        and print an optimal number of clusters     '''
    last = L[ - n_max_clusters :, 2 ]
    idxs = np.arange( 1, len(last) + 1)
    plt.plot( idxs, last[ : : - 1 ], label = 'distance from root for each cluster' )
    acceleration = np.diff( last, 2 )  # 2nd derivative of the distances
    acceleration_rev = acceleration[ : : -1 ]
    plt.plot( idxs[ : - 2 ] + 1, acceleration_rev, label = 'second derivative of "distance from root for each cluster"' )
    plt.show()
    k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
    plt.xlabel( 'number of cluster')
    plt.ylabel( 'distance from root' )
    plt.legend( )
    print( "optimal number of clusters:", k )


# In[ ]:


def Clustering( df, axis = 1, cluster_name = 'unnamedData_', clustering_type = 'optimal', is_relation_matrix = True, n_recursion_lim = 10, thres_skewed_ratio = 5, cluster_start = 0, min_num_genes = None, n_max_cluster = 10, linkage_method = None, cluster_method = 'maxclust', print_message = False ) :
    ''' Wrapper method of 'recursive_clustering' (recursively cluster and return subcluster-ordered entries) and 'hierarchical_clustering' method (return optimal ordering of entires). 
    set 'clustering_type' to 'recursive' or 'rc' for the former method and 'optimal' or 'hc' for the latter
    if list of df is given through 'df' argument, merge the dataframes (DataFrames has to have same indices though they can be in different order) and Clustering the DataFrame
    For asymmetric data (with different columns and indices), If 'axis' is 0, cluster columns. if 'axis' is 1, cluster indices. '''
    list_df = None # set a default value of list_df, which will contain list of DataFrame if it has been given through df argument
    if type( df ) is not pd.DataFrame :
        cluster_column_too = is_relation_matrix # store the value of 'is_relation_matrix' separately so that it can be used when returning the clustered result
        list_df, is_relation_matrix = deepcopy( df ), False # set 'is_relation_matrix' False if list of DataFrame has been given so that only indices are used for clustering
        df = list_df[ 0 ] # retrive the first dataframe in the given list of DataFrames 
        df = df.join( other = list_df[ 1 ], how = 'left', lsuffix = '_1', rsuffix = '_2' ) # join the first and second DataFrame, while adding suffix to column labels
        for index_df, a_df in zip( np.arange( 3, len( list_df ) + 1 ), list_df[ 2 : ] ) :
            df = df.join( other = a_df, how = 'left', rsuffix = '_' + str( index_df ) ) # join the thrid DataFrames and more DataFrames 
    if axis == 0 : # For asymmetric data (with different columns and indices), If 'axis' is 0, cluster columns. if 'axis' is 1, cluster indices.
        df = df.T
    if clustering_type == 'optimal' or clustering_type == 'hc' :
        dict_tree = hierarchical_clustering( df, cluster_name = cluster_name, is_relation_matrix = is_relation_matrix, cluster_start = cluster_start, min_num_genes = min_num_genes, n_max_cluster = n_max_cluster, linkage_method = linkage_method, cluster_method = cluster_method, print_message = print_message ) 
    else :
        if min_num_genes is None :
            min_num_genes = int( len( df ) * 0.01 ) if len( df ) > 2000 else 20 # set an optimal threshold for number of genes of a cluster to be subclustered (1% of all genes, minimum = 20 genes) 
        dict_tree = recursive_clustering( df, cluster_name = cluster_name, n_recursion_lim = n_recursion_lim, thres_skewed_ratio = thres_skewed_ratio, is_relation_matrix = is_relation_matrix, cluster_start = cluster_start, min_num_genes = min_num_genes, n_max_cluster = n_max_cluster, linkage_method = linkage_method, cluster_method = cluster_method, print_message = print_message ) 
    clustered_IDs = dict_tree[ 'Gene_IDs' ]
    if list_df is not None :
        if cluster_column_too :
            list_df_clustered = list( df.loc[ clustered_IDs, clustered_IDs ] for df in list_df )
        else :
            list_df_clustered = list( df.loc[ clustered_IDs ] for df in list_df )
        return list_df_clustered, dict_tree # return a list of Clustered dataFrames
    elif is_relation_matrix :
        return df.loc[ clustered_IDs, clustered_IDs ], dict_tree # return one Clustered relation_matrix dataFrame
    else :
        if axis == 0 : # For asymmetric data (with different columns and indices), If 'axis' is 0, cluster columns. if 'axis' is 1, cluster indices.
            return df.loc[ clustered_IDs ].T, dict_tree # return one Clustered dataFrame
        else :
            return df.loc[ clustered_IDs ], dict_tree # return one Clustered dataFrame


# In[ ]:


def CALCULATE_Correlation_Matrix( df, sample_list = None, Gene_Set = None, correlation_type = 'spearman', n_std_for_outliers = None, min_observation = 'auto' ) :
    ''' Calculate Spearman Correlation Coefficient for each gene-gene pairs and return return df_correl_mat
    default Gene_Set : indices of df (df.index.values)
    correlation_type = { 'pearson', 'kendall', 'spearman' } or lambda function that receive two numpy array
    it is expected to use require really long time computing when excluding outliers by using 'n_std_for_outliers' to calculate more robust correlation coefficients.  '''
    if Gene_Set is None : # set default Gene_Set
        Gene_Set = set( df.index.values )
    if sample_list is None : # set default sample_list
        sample_list = set( df.columns.values )
    valid_Gene_Set = set( df.index.values ).intersection( set( Gene_Set ) ) # valid Gene_Set is a set of genes that exist in df
    df_for_correl_matrix = df.loc[ valid_Gene_Set, sample_list ].dropna( ).transpose( ) # prepare a matrix for correlation matrix
    Gene_IDs = df_for_correl_matrix.columns.tolist( )
    data = df_for_correl_matrix.values
    if min_observation == 'auto' :
        min_observation = int( len( df_for_correl_matrix ) / 2 )
    if n_std_for_outliers is not None : # if 'n_std_for_outliers' has been given, create a mask for outliers for the use in subsequent calculation
        df_for_correl_matrix[ OUTLIERS_GET_mask_for_outliers( data.T, n_std_for_outliers = n_std_for_outliers ).T ] = np.nan # create a mask for outliers and replace outliers with np.nan 
    if correlation_type == 'spearman' : # if correlation_type == 's', compute Spearman's correl coefficients
        if n_std_for_outliers is not None :
            df_correl_mat = df_for_correl_matrix.corr( method = correlation_type, min_periods = min_observation )
        else :
            df_correl_mat = pd.DataFrame( stats.spearmanr( data ).correlation, index = Gene_IDs, columns = Gene_IDs ) # create DataFrame of correlation matrix 
    elif correlation_type == 'pearson' : # if correlation_type == 'p', compute Perason's correl coefficients
        if n_std_for_outliers is not None :
            df_correl_mat = df_for_correl_matrix.corr( method = correlation_type, min_periods = min_observation )
        else :
            df_correl_mat = pd.DataFrame( np.corrcoef( data.T ), index = Gene_IDs, columns = Gene_IDs ) # create DataFrame of correlation matrix 
    else :
        df_correl_mat = df_for_correl_matrix.corr( method = correlation_type, min_periods = min_observation ) # create DataFrame of correlation matrix 
    return df_correl_mat # return correlation matrix


# In[ ]:


def Correlation_Clustering( df, sample_list = None, Gene_Set = None, correlation_type = 'spearman', n_std_for_outliers = None, min_observation = 'auto', cluster_name = 'unnamedData_', is_relation_matrix = True, min_num_genes = 20, n_max_cluster = 10, clustering_type = 'optimal', linkage_method = None, cluster_method = 'maxclust', print_message = False ) :
    ''' Wrapper of two function : 'CALCULATE_Correlation_Matrix' and recursive function "Clustering". Caution has to be given that Clustering results might give wrongful insight that there are clear clusters, while there is actually no distinct clusters 
    df is a subset of df_proteome or df_phosphoproteome, or similar data frame with similar index and column structures
    return df_clustered_correl_mat, dict_tree_clusters as a tuple
    default Gene_Set : indices of df (df.index.values) '''
    df_correl_mat = CALCULATE_Correlation_Matrix( df = df, sample_list = sample_list, Gene_Set = Gene_Set, correlation_type = correlation_type, n_std_for_outliers = n_std_for_outliers, min_observation = min_observation ).dropna( ) # calculate correlation matrix
    _, dict_tree_clusters = Clustering( df_correl_mat, cluster_name = cluster_name, min_num_genes = min_num_genes, n_max_cluster = n_max_cluster, linkage_method = linkage_method, cluster_method = cluster_method, print_message = print_message, clustering_type = clustering_type ) # perform a recursive clustering of df
    df_clustered_correl_mat = df_correl_mat.loc[ dict_tree_clusters[ 'Gene_IDs' ], dict_tree_clusters[ 'Gene_IDs' ] ] # rearrange indices and columns of df to cluster genes according to the result
    return df_clustered_correl_mat, dict_tree_clusters


# In[ ]:


def ORANGE_CLUSTERING_explore_tree_node( tree_node, list_left_to_right ) :
    ''' a recursive method that travels Orange3's cluster node from left to right and retrive the order of index in the tree '''
    if tree_node.is_leaf :
        list_left_to_right.append( tree_node.value.index )
        return list_left_to_right
    else :
        list_left_to_right = ORANGE_CLUSTERING_explore_tree_node( tree_node.left, list_left_to_right )
        list_left_to_right = ORANGE_CLUSTERING_explore_tree_node( tree_node.right, list_left_to_right )
        return list_left_to_right


# In[ ]:


"""def hierarchical_clustering( df = None, cluster_name = 'unnamedData_' , is_relation_matrix = True, cluster_start = 0, min_num_genes = 20, n_max_cluster = 10, linkage_method = None, cluster_method = 'maxclust', print_message = True ):
    '''    Optimal leaf ordering by using Orange3 library    '''
    Gene_IDs_before_clustering = list( df.index.values ) # retrive indices (Gene_IDs)'
    dist_matrix = od.Euclidean( df.values ) # calculate distance matrix
    L = och.dist_matrix_linkage( dist_matrix ) # calculate 
    och.cophenetic_correlation(cluster, matrix)
    tree_root_node = och.tree_from_linkage( L )
    ordered_tree = och.optimal_leaf_ordering( tree_root_node, dist_matrix )
    list_optimal_ordering = ORANGE_CLUSTERING_explore_tree_node( ordered_tree, list( ) )
    df_clustered = df.iloc[ list_optimal_ordering, list_optimal_ordering ]
    return dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = subclustered_Gene_IDs, subcluster_tree = dict_subclusters, Clustering_Cophenetic_Coeff = Clustering_Cophenetic_Coeff, linkage_method = linkage_method, cluster_method = cluster_method )"""


# In[ ]:


# def hierarchical_clustering( df = None, cluster_name = 'unnamedData_' , is_relation_matrix = True, cluster_start = 0, min_num_genes = 20, n_max_cluster = 10, linkage_method = None, cluster_method = 'maxclust', print_message = True ):
#     '''    Optimal leaf ordering by using Orange3 library. use scipy module to optimize linkage method if linkage_method = None. See 'recursive_clustering' method for details. 
#     'cluster_method', 'n_max_cluster', 'cluster_start', 'min_num_genes' arguments are not functioning (they exist to increase code compatibility with recursive_clustering method)   '''
#     Gene_IDs_before_clustering = df.index.values # retrive indices (Gene_IDs)'
#     dist_matrix = od.Euclidean( df.values ) # calculate distance matrix
#     condensed_d = och.condensedform( dist_matrix ) # retrive condensed form of distance matrix, which will be used for scipy Cophenetic_Coeffs computation
#     if print_message :
#         print( cluster_name, '(Step 1) distances calculated ({n_genes} Genes)'.format( n_genes = len( Gene_IDs_before_clustering ) ) ) # report progress of clustering step (2)
#     if linkage_method is None : # if linkage_method has not been specified, try all methods and use method that gives optimal Clustering_Cophenetic_Coeff
#         linkage_methods = [ 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' ] # define all available linkage_methods
#         list_Clustering_Cophenetic_Coeffs = list( )
#         dict_method_L = dict( )
#         for method in linkage_methods : # for each method, calculate and store Linkage map and Cophenetic_Coeff of Clustering
#             dict_method_L[ method ] = sch.linkage( condensed_d, method = method ) # clustering using a specified method 
#             list_Clustering_Cophenetic_Coeffs.append( sch.cophenet( dict_method_L[ method ], condensed_d )[ 0 ] ) # add calculated Clustering_Cophenetic_Coeff to the array
#         optimal_method_index = np.argmax( list_Clustering_Cophenetic_Coeffs ) # find index of optimal linkage_method by finding positon of maximum Clustering_Cophenetic_Coeff (closest to 1)
#         linkage_method = linkage_methods[ optimal_method_index ] # retrive optimal linkage_method, linkage_map, and Clustering_Cophenetic_Coeff
#         Clustering_Cophenetic_Coeff = list_Clustering_Cophenetic_Coeffs[ optimal_method_index ]
#         L = dict_method_L[ linkage_method ]
#         if print_message :
#             print( cluster_name, '(Step 2) optimal linkage_map calculation method found :', linkage_method ) # report progress of clustering step (2)
#     tree_root_node = och.tree_from_linkage( L )
#     ordered_tree = och.optimal_leaf_ordering( tree_root_node, dist_matrix )
#     list_optimal_ordering = ORANGE_CLUSTERING_explore_tree_node( ordered_tree, list( ) )
#     Gene_IDs_after_clustering = Gene_IDs_before_clustering[ list_optimal_ordering ]
#     return dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = Gene_IDs_after_clustering, subcluster_tree = None, Clustering_Cophenetic_Coeff = Clustering_Cophenetic_Coeff, linkage_method = linkage_method, cluster_method = None )


# In[ ]:


# def recursive_clustering( df = None, df_d = None, cluster_name = 'unnamedData_', n_recursion_lim = 10, thres_skewed_ratio = 5, is_relation_matrix = True, cluster_start = 0, min_num_genes = 20, n_max_cluster = 10, linkage_method = None, cluster_method = 'maxclust', print_message = True ):
#     ''' Recursively clustering correlation matrix until number of gene in a subcluster is less than a given minimal number of genes
#     df = Dataframe of correlation matrix or any matrix of which indices and columns are same 
#     min_num_genes = 100
#     is_relation_matrix : if True, assume the given 'df' is a DataFame of which columns and indices are Gene_IDs (the same number of cols and rows). if False, assume Gene_IDs are indices (* can reduce computation time for relation matrices)
#     cluster_method : 'distance' or 'maxclust' (elbow method)
#                     if 'maxclust' method (elbow method) yield only one actual cluster, change cluster_method to  'distance' method
#     linkage_method : mostly 'centroid' or 'complete' works best for correlation matrix, and 'single' works worst
#                      By default (when None is given), automtically choose likage_method that yield optimal Clustering_Cophenetic_Coeff (that most close to 1)
#     cluster_name convention : 4.1.2.3 (subcluster numbers seperated by columns)
#     return a tree structure of sub clusters in each clusters as a dictionary :
#         node structure : dict( cluster_name = , cluster_start = , Gene_IDs = , subcluster_tree = , Clustering_Cophenetic_Coeff = , linkage_method = , cluster_method = ) : subcluster_tree is None if the node is the leaf node
#     Clustering_Cophenetic_Coeff : Cophenetic Correlation Coefficient of clustering; the closer the value is to 1, the better the clustering preserves the original distances
#     pre-computed pairwise distance can be given through 'df_d' argument   '''
#     FLAG_auto_linkage_method = False # initial FLAG for setting subcluster linkage method
#     Gene_IDs_before_clustering, n_genes = list( df.index.values ), len( df.index.values ) # retrive indices (Gene_IDs)'
#     dict_leaf_cluster = dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = Gene_IDs_before_clustering, subcluster_tree = None ) # set a leaf node of cluster tree
#     if len( Gene_IDs_before_clustering ) < min_num_genes or n_recursion_lim < cluster_name.count( '.' ): # a terminating condition : return a leaf node when number of genes is smaller than a given threshold or number of recursion exceeded the recursion limit
#         return dict_leaf_cluster
#     if df_d is None or len( df_d ) != len( df ) : # if distance matrix was not given or its dimension is incompatible with a given DataFrame, calculate distance matrix
#         d = sch.distance.pdist( df.values ) # condensed form of pairwise distance matrix (Euclidean)
#         df_d = pd.DataFrame( sch.distance.squareform( d ), index = Gene_IDs_before_clustering, columns = Gene_IDs_before_clustering ) # retrive the square form of distance matrix
#     else :
#         d = sch.distance.squareform( df_d.values ) # retrive the condensed form of distance matrix
#     if print_message :
#         print( cluster_name, '(Step 1) distances calculated ({} Genes)'.format( n_genes ) ) # report progress of clustering step (2)
#     if linkage_method is None : # if linkage_method has not been specified, try all methods and use method that gives optimal Clustering_Cophenetic_Coeff
#         FLAG_auto_linkage_method = True
#         linkage_methods = [ 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' ] # define all available linkage_methods
#         list_Clustering_Cophenetic_Coeffs = list( )
#         dict_method_L = dict( )
#         for method in linkage_methods : # for each method, calculate and store Linkage map and Cophenetic_Coeff of Clustering
#             dict_method_L[ method ] = sch.linkage( d, method = method ) # clustering using a specified method 
#             list_Clustering_Cophenetic_Coeffs.append( sch.cophenet( dict_method_L[ method ], d )[ 0 ] ) # add calculated Clustering_Cophenetic_Coeff to the array
#         optimal_method_index = np.argmax( list_Clustering_Cophenetic_Coeffs ) # find index of optimal linkage_method by finding positon of maximum Clustering_Cophenetic_Coeff (closest to 1)
#         linkage_method = linkage_methods[ optimal_method_index ] # retrive optimal linkage_method, linkage_map, and Clustering_Cophenetic_Coeff
#         Clustering_Cophenetic_Coeff = list_Clustering_Cophenetic_Coeffs[ optimal_method_index ]
#         L = dict_method_L[ linkage_method ]
#         del dict_method_L # delete Linkage data to free memory
#         if print_message :
#             print( cluster_name, '(Step 2) optimal linkage_map calculation method found :', linkage_method ) # report progress of clustering step (2)
#     else : # if linkage_method has been specified, use the method to calculate linkage_map and Clustering_Cophenetic_Coeff
#         L = sch.linkage( d, method = linkage_method )
#         Clustering_Cophenetic_Coeff = sch.cophenet( L, d )[ 0 ] # Cophenetic Correlation Coefficient of clustering 
#     if cluster_method == 'maxclust' or cluster_method == 'elbow' : # define optimal number of clusters using elbow method
#         acceleration_rev = np.diff( L[ -n_max_cluster:, 2 ], 2 )[ : : -1 ]  # reverse of 2nd derivative of the last n_max_cluster number of distances
#         opt_n_clusters = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters # set optimal number of clusters
#         ind = sch.fcluster( L, opt_n_clusters, criterion = 'maxclust' )
#         n_clusters = ind.max( ) # retrive actual number of clusters
#         skewed_ratio = ( ( LIST_COUNT( ind, duplicate_filter = None ).values / n_genes ) ** 2 ).sum( ) * n_clusters # calculate skewed ratio, which is close to 1 if a cluster was equally splitted subclusters and becomes larger as one subcluster take most genes of the cluster
#         if print_message : 
#             print( cluster_name, '(Step 2.5) elbow method performed, and optimal number of clusters : {}, actual number of clusters : {}, skewed ratio : {:.3f}'.format( opt_n_clusters, n_clusters, skewed_ratio ) )
#         if n_clusters == 1 or thres_skewed_ratio < skewed_ratio : # if 'maxclust' method (elbow method) yield only one actual cluster, change cluster_method to  'distance' method
#             cluster_method = 'distance'
#             print( cluster_name, "one cluster possible by elbow method or skewed ratio is higher than the threshold, change clustering method to 'distance'" )
#     if cluster_method == 'distance' : # as alternatively, use maximum distance to set number of clusters
#         d_max = d.max( )
#         d_max_ratio = 0.5
#         ind = sch.fcluster( L, d_max_ratio * d_max, criterion = 'distance' )
#         n_clusters = ind.max( ) # retrive the actual number of clusters
#         skewed_ratio = ( ( LIST_COUNT( ind, duplicate_filter = None ).values / n_genes ) ** 2 ).sum( ) * n_clusters
#         while thres_skewed_ratio < skewed_ratio :
#             d_max_ratio = d_max_ratio * 0.80 # reduce maximum distance required for clustering by 20% in each iteration
#             ind = sch.fcluster( L, d_max_ratio * d_max, criterion = 'distance' ) # cluster with new maximum distance and re-calculate skewed_ratio
#             n_clusters = ind.max( ) # retrive the actual number of clusters
#             skewed_ratio = ( ( LIST_COUNT( ind, duplicate_filter = None ).values / n_genes ) ** 2 ).sum( ) * n_clusters
#             if d_max_ratio < 1e-4 or n_clusters * 20 > n_genes : # if distance become too small or n_clusters becomes to large, stop clustering 
#                 break
#         if n_clusters == 1 : # SAFE TERMINATING CONDITION : if there is only one cluster possible, since it satisfies the terminating conditions, return leaf_cluster_node
#             print( cluster_name, 'one cluster possible, and clustering stopped')
#             return dict_leaf_cluster 
#     if print_message : 
#         print( cluster_name, '(Step 3) clustering performed, number of clusters : {}, Clustering_Cophenetic_Coeff : {:.3f}, skewed ratio : {:.3f}'.format( n_clusters, Clustering_Cophenetic_Coeff, skewed_ratio ) ) # print Clustering_Cophenetic_Coeff
#     del d # delete distance vector to free its memory
#     del L # delete linkage map to free its memory
#     clustered_Gene_IDs = [ Gene_IDs_before_clustering[ i ] for i in list( np.argsort( ind ) ) ] # sort column names according to the clustering result
#     df = df.loc[ clustered_Gene_IDs, clustered_Gene_IDs ] if is_relation_matrix else df.loc[ clustered_Gene_IDs ] # re-aligning columns (Gene_IDs) and/or indices according to the clustering result
#     df_d = df_d.loc[ clustered_Gene_IDs, clustered_Gene_IDs ] # re-aligning columns (Gene_IDs) and indices of distance matrix according to the clustering result
#     dict_subclusters = dict( ) # create a dictionary for subclusters
#     # recursive clustering for each newly assigned cluster
#     sub_clus_linkage_method = None if FLAG_auto_linkage_method else linkage_method # set subcluster linkage method    
#     subcluster_sizes = pd.Series( dict( collections.Counter( ind ) ) ).sort_index( ).values # retrive list of cluster sizes in the same order as cluster assignment number
#     del ind # delete clustered result to free its memory
#     subcluster_numbers = [ i + 1 for i in range( len( subcluster_sizes ) ) ] # list of assigned cluster numbers
#     subcluster_start = 0 # set the position of the first cluster start site as 0 (start of an array)
#     subclustered_Gene_IDs = [ ] # list of Gene_IDs that were subclustered recursively
#     for subcluster_size, subcluster_number in zip( subcluster_sizes, subcluster_numbers ) :
#         subcluster_Gene_IDs = clustered_Gene_IDs[ subcluster_start : subcluster_start + subcluster_size ] # retrive columns (Gene_IDs) of a sub-cluster
#         df_subcluster = df.loc[ subcluster_Gene_IDs, subcluster_Gene_IDs ] if is_relation_matrix else df.loc[ subcluster_Gene_IDs ] # create a dataframe for a sub-cluster
#         df_d_subcluster = df_d.loc[ subcluster_Gene_IDs, subcluster_Gene_IDs ] # create a dataframe of pairwise distance for a sub-cluster
#         subcluster_name = cluster_name + '.' + str( subcluster_number ) # set a sub-cluster name 
#         # recursively cluster genes of a subcluster
#         dict_subclusters[ subcluster_name ] = recursive_clustering( df_subcluster, df_d = df_d_subcluster, cluster_name = subcluster_name, n_recursion_lim = n_recursion_lim, is_relation_matrix = is_relation_matrix, cluster_start = cluster_start + subcluster_start, min_num_genes = min_num_genes, n_max_cluster = n_max_cluster, linkage_method = sub_clus_linkage_method, cluster_method = cluster_method, print_message = print_message ) 
#         subclustered_Gene_IDs.extend( dict_subclusters[ subcluster_name ][ 'Gene_IDs' ] ) # extend a list of clustered gene_ids of a subcluster to the subclustered gene_id list
#         subcluster_start += subcluster_size # set the next sub-cluster start position
#     return dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = subclustered_Gene_IDs, subcluster_tree = dict_subclusters, Clustering_Cophenetic_Coeff = Clustering_Cophenetic_Coeff, linkage_method = linkage_method, cluster_method = cluster_method )


# In[ ]:


""" # frozen at 20190602
def recursive_clustering( df = None, cluster_name = 'unnamedData_' , is_relation_matrix = True, cluster_start = 0, min_num_genes = 20, n_max_cluster = 10, linkage_method = None, cluster_method = 'maxclust', print_message = True ):
    '''
    Recursively clustering correlation matrix until number of gene in a subcluster is less than a given minimal number of genes
    df = Dataframe of correlation matrix or any matrix of which indices and columns are same 
    min_num_genes = 100
    is_relation_matrix : if True, assume the given 'df' is a DataFame of which columns and indices are Gene_IDs (the same number of cols and rows). if False, assume Gene_IDs are indices (* can reduce computation time for relation matrices)
    cluster_method : 'distance' or 'maxclust' (elbow method)
                    if 'maxclust' method (elbow method) yield only one actual cluster, change cluster_method to  'distance' method
    linkage_method : mostly 'centroid' or 'complete' works best for correlation matrix, and 'single' works worst
                     By default (when None is given), automtically choose likage_method that yield optimal Clustering_Cophenetic_Coeff (that most close to 1)
    cluster_name convention : 4.1.2.3 (subcluster numbers seperated by columns)
    return a tree structure of sub clusters in each clusters as a dictionary :
        node structure : dict( cluster_name = , cluster_start = , Gene_IDs = , subcluster_tree = , Clustering_Cophenetic_Coeff = , linkage_method = , cluster_method = ) : subcluster_tree is None if the node is the leaf node
    Clustering_Cophenetic_Coeff : Cophenetic Correlation Coefficient of clustering; the closer the value is to 1, the better the clustering preserves the original distances
    '''
    Gene_IDs_before_clustering = list( df.index.values ) # retrive indices (Gene_IDs)'
    dict_leaf_cluster = dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = Gene_IDs_before_clustering, subcluster_tree = None ) # set a leaf node of cluster tree
    if len( Gene_IDs_before_clustering ) < min_num_genes : # a terminating condition : return a leaf node when number of genes is smaller than a given threshold
        return dict_leaf_cluster
    d = sch.distance.pdist( df.values )   # vector of pairwise distances
    if print_message :
        print( cluster_name, '(Step 1) distances calculated ({n_genes} Genes)'.format( n_genes = len( Gene_IDs_before_clustering ) ) ) # report progress of clustering step (2)
    if linkage_method is None : # if linkage_method has not been specified, try all methods and use method that gives optimal Clustering_Cophenetic_Coeff
        linkage_methods = [ 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' ] # define all available linkage_methods
        list_Clustering_Cophenetic_Coeffs = list( )
        dict_method_L = dict( )
        for method in linkage_methods : # for each method, calculate and store Linkage map and Cophenetic_Coeff of Clustering
            dict_method_L[ method ] = sch.linkage( d, method = method ) # clustering using a specified method 
            list_Clustering_Cophenetic_Coeffs.append( sch.cophenet( dict_method_L[ method ], d )[ 0 ] ) # add calculated Clustering_Cophenetic_Coeff to the array
        optimal_method_index = np.argmax( list_Clustering_Cophenetic_Coeffs ) # find index of optimal linkage_method by finding positon of maximum Clustering_Cophenetic_Coeff (closest to 1)
        linkage_method = linkage_methods[ optimal_method_index ] # retrive optimal linkage_method, linkage_map, and Clustering_Cophenetic_Coeff
        Clustering_Cophenetic_Coeff = list_Clustering_Cophenetic_Coeffs[ optimal_method_index ]
        L = dict_method_L[ linkage_method ]
        if print_message :
            print( cluster_name, '(Step 2) optimal linkage_map calculation method found :', linkage_method ) # report progress of clustering step (2)
    else : # if linkage_method has been specified, use the method to calculate linkage_map and Clustering_Cophenetic_Coeff
        L = sch.linkage( d, method = linkage_method )
        Clustering_Cophenetic_Coeff = sch.cophenet( L, d )[ 0 ] # Cophenetic Correlation Coefficient of clustering 
    if cluster_method == 'maxclust' or cluster_method == 'elbow' : # define optimal number of clusters using elbow method
        acceleration_rev = np.diff( L[ -n_max_cluster:, 2 ], 2 )[ : : -1 ]  # reverse of 2nd derivative of the last n_max_cluster number of distances
        opt_n_clusters = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters # set optimal number of clusters
        if print_message : 
            print( cluster_name, '(Step 2.5) elbow method performed, and optimal number of clusters :', opt_n_clusters )
        ind = sch.fcluster( L, opt_n_clusters, criterion = 'maxclust' )
        if ind.max( ) == 1 : # if 'maxclust' method (elbow method) yield only one actual cluster, change cluster_method to  'distance' method
            cluster_method = 'distance'
            print( cluster_name, "one cluster possible by elbow method, change clustering method to 'distance'" )
    if cluster_method == 'distance' : # as alternatively, use maximum distance to set number of clusters
        ind = sch.fcluster( L, 0.5 * d.max( ), criterion = 'distance' )
    k = ind.max( ) # retrive the actual number of clusters
    if k < 2 : # SAFE TERMINATING CONDITION : if there is only one cluster possible, since it satisfies the terminating conditions, return leaf_cluster_node
        print( cluster_name, 'one cluster possible, and clustering stopped')
        return dict_leaf_cluster 
    if print_message : 
        print( cluster_name, '(Step 3) clustering performed, number of clusters :', k, '\t Clustering_Cophenetic_Coeff :', Clustering_Cophenetic_Coeff ) # print Clustering_Cophenetic_Coeff
    clustered_Gene_IDs = [ Gene_IDs_before_clustering[ i ] for i in list( np.argsort( ind ) ) ] # sort column names according to the clustering result
    if is_relation_matrix :
        df = df.loc[ clustered_Gene_IDs, clustered_Gene_IDs ] # re-aligning columns (Gene_IDs) and indices according to the clustering result
    else :
        df = df.loc[ clustered_Gene_IDs ] # re-aligning only indices (Gene_IDs) according to the clustering result if it is not an relational matrix
    dict_subclusters = dict( ) # create a dictionary for subclusters
    # recursive clustering for each newly assigned cluster
    subcluster_sizes = pd.Series( dict( collections.Counter( ind ) ) ).sort_index( ).values # retrive list of cluster sizes in the same order as cluster assignment number
    subcluster_numbers = [ i + 1 for i in range( len( subcluster_sizes ) ) ] # list of assigned cluster numbers
    subcluster_start = 0 # set the position of the first cluster start site as 0 (start of an array)
    subclustered_Gene_IDs = [ ] # list of Gene_IDs that were subclustered recursively
    for subcluster_size, subcluster_number in zip( subcluster_sizes, subcluster_numbers ) :
        subcluster_Gene_IDs = clustered_Gene_IDs[ subcluster_start : subcluster_start + subcluster_size ] # retrive columns (Gene_IDs) of a sub-cluster
        if is_relation_matrix :
            df_subcluster = df.loc[ subcluster_Gene_IDs, subcluster_Gene_IDs ] # create a dataframe for a sub-cluster of relation matrix (subsetting both row and column thus reducing data points can increase computation time)
        else :
            df_subcluster = df.loc[ subcluster_Gene_IDs ] # create a dataframe for a sub-cluster
        subcluster_name = cluster_name + '.' + str( subcluster_number ) # set a sub-cluster name 
        # recursively cluster genes of a subcluster
        dict_subclusters[ subcluster_name ] = recursive_clustering( df_subcluster, cluster_name = subcluster_name, is_relation_matrix = is_relation_matrix, cluster_start = cluster_start + subcluster_start, min_num_genes = min_num_genes, n_max_cluster = n_max_cluster, linkage_method = linkage_method, cluster_method = cluster_method, print_message = print_message ) 
        subclustered_Gene_IDs.extend( dict_subclusters[ subcluster_name ][ 'Gene_IDs' ] ) # extend a list of clustered gene_ids of a subcluster to the subclustered gene_id list
        subcluster_start += subcluster_size # set the next sub-cluster start position
    return dict( cluster_name = cluster_name, cluster_start = cluster_start, Gene_IDs = subclustered_Gene_IDs, subcluster_tree = dict_subclusters, Clustering_Cophenetic_Coeff = Clustering_Cophenetic_Coeff, linkage_method = linkage_method, cluster_method = cluster_method )"""


# In[ ]:


def print_cluster_tree( cluster_tree, print_cluster_location = True, threshold = 50, d = 0 ):
    '''
    function to visualize recursive clustering result 
    recursive method to visualize a dictionary-based tree
    threshold = 50 : threshold of number of genes for print information about current cluster
    '''
    if len( cluster_tree[ 'Gene_IDs'] ) < threshold : # if the number of gene of current cluster is less then a threshold, skip printing
        return
    if print_cluster_location : # print information about current cluster 
        print( "\t" * d, cluster_tree[ 'cluster_name' ], 'number of genes :', len( cluster_tree[ 'Gene_IDs'] ), '\tstart and end position =', cluster_tree[ 'cluster_start' ], ':', cluster_tree[ 'cluster_start' ] + len( cluster_tree[ 'Gene_IDs' ] ) )
    else :
        print( "\t" * d, cluster_tree[ 'cluster_name' ], 'number of genes :', len( cluster_tree[ 'Gene_IDs'] ) )
    if cluster_tree[ 'subcluster_tree' ] is not None : # recursively plot subclusters
        for key, val in cluster_tree[ 'subcluster_tree' ].items():
            print_cluster_tree( val, print_cluster_location = print_cluster_location, threshold = threshold, d = d + 1 ) # recursion


# In[ ]:


bbox_props_white = dict( boxstyle = "round", fc = "w", ec="0.5", alpha = 0.7 )  # define box properties that will be used in a function below for annotation
def Plot_corr_matrix__Gene_IDs( df, df_data = None, color_limit = 'auto', color_limit_infobar = None, color_map_diversing = True, show_inforbar_separately = True, color_map = 'bwr', Gene_IDs_interest = None, mask_thres = None, sample_list_min_by_max = None, add_Log2FC = False, n_factors = None, save_fig = False, show_grid = False, show_colorbar = True, Genes_of_interest_name = 'unknown_gene_set', show_title = True, mask_alpha = 0.4, image_dpi = 200, df_name = 'Data', Log2FC_name = 'Log2FC\nTumor to Normal\n', gap_width = 1, title = '', folder = None, gene_list_annotation = None, gene_annotation_text_size = 10 ) :
    ''' Plot a graphical correlation matrix* for a correlation matrix dataframe. 
        * a correlation matrix means index and columns have same labels, and element at i, j indicates some sort of relationship between i and j)
    locate the position of a given list of gene_ids by crosses of black lines,
    or locate the positions where value is higher than the given threshold.
    Also, visualize Log2FoldChange values along with correlation matrix
    visuallze matrix tightly without gene names to speed up the plotting.
    mask_alpha = 0.5 : alpha value for the mask for locating genes
    df : dataframe or a list of DataFrames
    df_data : dataframe from which Log2FC will be calculated
    df_name : string or list of strings of df_names
    color_map : a name of colormap or list of names of color_maps for plots 
    color_limit : value for upper and lower limit for colormap values. if 'color_limit' = 'auto', set color_limit as mean + std * 1.2
    gene_list_annotation : see Function 'GENE_Annotate_df_for_plot'    '''
    if ( ( add_Log2FC or n_factors is not None ) or sample_list_min_by_max is not None ) and df_data is None : # set a defualt dataframe to create infobar and use them for ploting 
        df_data = df_proteome_unshared
    Symbols_annotation, list_positions_annotation = GENE_Annotate_df_for_plot( df, gene_list_annotation ) # retrive symbols and positions of genes that will be annotated on the plot
    if type( df ) is list :
        if len( df ) != len( df_name ) : # checking and prepare drawing multiple plots 
            print( 'number of df and df_name does not match' )
            return - 1
        if type( color_map ) is not list :
            color_map = [ color_map ] * len( df )
        elif len( color_map ) != len( df ) :
            print( 'number of df and colormap does not match' )
            return - 1
        fig, ax = plt.subplots( 1, len( df ), sharex = True, sharey = True )
        for a_ax, a_df, a_df_name, a_color_map in zip( ax, df, df_name, color_map ) :
            a_ax = subplot_relation_matrix__Gene_IDs( fig, a_ax, a_df, df_data = df_data, color_limit = color_limit, color_limit_infobar = color_limit_infobar, color_map_diversing = color_map_diversing, color_map = a_color_map, Gene_IDs_interest = Gene_IDs_interest, mask_thres = mask_thres, sample_list_min_by_max = sample_list_min_by_max, add_Log2FC = add_Log2FC, n_factors = n_factors, show_grid = show_grid, show_colorbar = show_colorbar, Genes_of_interest_name = Genes_of_interest_name, show_title = show_title, mask_alpha = mask_alpha, df_name = a_df_name, Log2FC_name = Log2FC_name, gap_width = gap_width )
    else :
        if show_inforbar_separately and add_Log2FC : # if 'show_inforbar_separately' is true, create two axes and plot matrix and inforbar separately
            fig, axes = plt.subplots( 1, 2, sharex = True, sharey = True )
            ax_matrix, ax_infobar = axes
            ax = subplot_relation_matrix__Gene_IDs( fig, ax_matrix, df, ax_infobar = ax_infobar, df_data = df_data, color_limit = color_limit, color_limit_infobar = color_limit_infobar, color_map_diversing = color_map_diversing, color_map = color_map, Gene_IDs_interest = Gene_IDs_interest, mask_thres = mask_thres, sample_list_min_by_max = sample_list_min_by_max, add_Log2FC = add_Log2FC, n_factors = n_factors, show_grid = show_grid, show_colorbar = show_colorbar, Genes_of_interest_name = Genes_of_interest_name, show_title = show_title, mask_alpha = mask_alpha, df_name = df_name, Log2FC_name = Log2FC_name, gap_width = gap_width )
            plt.subplots_adjust( wspace = 0 ) # remove s white space between the two axes
        else :
            fig, ax = plt.subplots( 1, 1 )
            ax = subplot_relation_matrix__Gene_IDs( fig, ax, df, df_data = df_data, color_limit = color_limit, color_limit_infobar = color_limit_infobar, color_map_diversing = color_map_diversing, color_map = color_map, Gene_IDs_interest = Gene_IDs_interest, mask_thres = mask_thres, sample_list_min_by_max = sample_list_min_by_max, add_Log2FC = add_Log2FC, n_factors = n_factors, show_grid = show_grid, show_colorbar = show_colorbar, Genes_of_interest_name = Genes_of_interest_name, show_title = show_title, mask_alpha = mask_alpha, df_name = df_name, Log2FC_name = Log2FC_name, gap_width = gap_width )
        if gene_list_annotation is not None : # if 'gene_list_annotation' is given, annotated each Gene_Symbol on the plot 
            for Symbol, position in zip( Symbols_annotation, list_positions_annotation ) :
                ax.text( position, position, Symbol, ha = "right", va = "center", size = gene_annotation_text_size, bbox = bbox_props_white )
    if save_fig : # save and close figure if save_fig = True
        plt.savefig( folder + To_window_path_compatible_str( title ) + '.png', dpi = image_dpi )
        plt.close( )
    else :
        return fig, ax


# In[ ]:


def subplot_relation_matrix__Gene_IDs( fig, ax, df, ax_infobar = None, df_data = None, color_limit = 'auto', color_limit_infobar = None, color_map_diversing = True, color_map = 'bwr', Gene_IDs_interest = None, mask_thres = None, sample_list_min_by_max = None, add_Log2FC = True, n_factors = 3, show_grid = False, show_colorbar = True, Genes_of_interest_name = 'unknown_gene_set', show_title = True, mask_alpha = 0.4, df_name = 'df', Log2FC_name = 'Log2FC\nTumor to Normal\n', gap_width = 1 ) :
    ''' receive a figure and a subplot and plot a graphical correlation matrix* for a correlation matrix dataframe. Plot the matrix and infobar separately. Infobar (Log2FC or Factor Analysis) is plotted with 'bwr' colormap
    * a correlation matrix means index and columns have same labels, and element at i, j indicates some sort of relationship between i and j)
    locate the position of a given list of gene_ids by crosses of black lines,
    or locate the positions where value is higher than the given threshold.
    Also, visualize Log2FoldChange values along with correlation matrix
    visuallze matrix tightly without gene names to speed up the plotting.
    mask_alpha = 0.5 : alpha value for the mask for locating genes
    df : dataframe or list of DataFrame
    df_data : dataframe from which Log2FC will be calculated
    df_name : string or list of strings of df_names
    color_limit : value for upper and lower limit for colormap values. if 'color_limit' = 'auto', set color_limit as mean + std * 1.2
    n_factors : if an interger number is given, Perfoem Factor Analysis with a given number of factors and plot them along with log2FC
    gap_width : a width of white pixels separating Log2FC and Factor Analysis result from each other 
    color_limit_infobar : by default, color_limit_infobar = color_limit
    '''
    title = df_name # set initial title of the subplot
    corr = df.values # retrive correlation matrix
    n = len( corr ) # retrive number of entries
    if Gene_IDs_interest is not None : # if 'Gene_IDs_interest' is given, plot matrix and infobar together (plotting separately is currently not supported) 
        ax_infobar = None 
    if color_limit == 'auto' : #  if 'color_limit' = 'auto', automatically set color_limit as mean + std * 1.2
        abs_data = np.abs( df.values )
        color_limit = abs_data.mean( ) + abs_data.std( ) * 1.2
    if color_limit_infobar is None : # set color_limit as color_limit_infobar if its value has not been given
        color_limit_infobar = color_limit
    if add_Log2FC and df_data is None : # set defualt dataframe for create Log2FC and use them for ploting 
        df_data = df_proteome_unshared
    if ax_infobar is not None : # set values when inforbar is plotted separately
        arr_infobar = 'Empty Infobar'
        label_infobar = ''
        length_Log2FC_infobar = n
        n_factors = None # currently n_factors are not supported in two axes-subplot
    else :
        length_Log2FC_infobar = int( len( df ) / 10 )
        if n_factors is not None :
            dislay_length_for_a_factor = int( len( df ) / 20 )
            length_Factor_Analysis_infobar = dislay_length_for_a_factor * n_factors
    if add_Log2FC : # add Log2FC inforbar if 'add_Log2FC' is set to True
        arr_Log2FC = np.zeros( ( n, length_Log2FC_infobar ) ).T
        arr_mask = deepcopy( arr_Log2FC ) # an empty array that will be added to mask
        if sample_list_min_by_max is not None : # 'sample_list_min_by_max' is not None, override Log2FC infobar with min_by_max values
            arr_Log2FC[ gap_width : , : ] = Calculate_min_by_max( df_data, sample_list = sample_list_min_by_max, Gene_Set = df.index.values ).values # Calculate and broadcast min_by_max to the 2-D array for better visualization 
        else :
            arr_Log2FC[ gap_width : , : ] = Calculate_Log2FC_p_value__A_vs_B( df_data, Gene_IDs = df.index.values ).loc[ df.index.values, 'Log2_Fold_Change' ].values # Calculate and broadcast Log2FC to the 2-D array for better visualization 
        if ax_infobar is None :
            corr = np.vstack( ( corr.T, arr_Log2FC ) ).T # add Log2FC to correlation matrix
            ax.annotate( Log2FC_name, xy = ( n, 0 ), fontsize = 10 )
        else : # when inforbar is plotted separately
            arr_infobar = arr_Log2FC.T
            label_infobar += Log2FC_name
    if n_factors is not None : # add Factor Analysis inforbar if 'n_factors' is not None
        arr_FA = np.zeros( ( n, length_Factor_Analysis_infobar ) ).T
        arr_mask_FA = deepcopy( arr_FA ) # an empty array that will be added to mask
        df_FA = CALCULATE_Factor_Analysis__df( df_data, N_samples, n_factors = 5, Gene_Set = df.index.values ).loc[ df.index.values ] # Perform Factor Analysis and rearrange the result according to the given relation_matrix
        for index_factor, data_of_a_factor in zip( np.arange( n_factors ), df_FA.values.T ) :
            arr_FA[ gap_width + index_factor * dislay_length_for_a_factor : ( index_factor + 1 ) * dislay_length_for_a_factor, : ] = data_of_a_factor # broadcast data of a factor to the array that will be displayed
        label_Factor_Analysis = 'Factor Analysis ({n_factors} Factors)\n'.format( n_factors = n_factors )
        if ax_infobar is None :
            corr = np.vstack( ( corr.T, arr_FA ) ).T # add FA result to correlation matrix and Log2FC (if added)
            ax.annotate( label_Factor_Analysis, xy = ( n + int( dislay_length_for_a_factor * n_factors / 2 ), 0 ), fontsize = 10 )  # annotate FA result     
        else : # when inforbar is plotted separately
            label_infobar += label_Factor_Analysis
            if type( arr_infobar ) == 'Empty Infobar' :
                arr_infobar = arr_FA.T
            else :
                arr_infobar = np.vstack( ( arr_infobar.T, arr_FA ) ).T
    if ax_infobar is not None or ( not add_Log2FC and n_factors is None ) : # Plot the correlation matrix
        im = ax.imshow( corr, cmap = color_map )
    else :
        length_matrix, length_matrix_plus_infobar = np.shape( corr )
        mask_infobar = np.zeros_like( corr ).astype( bool )
        mask_infobar[ :, length_matrix : ] = True
        corr_matrix = np.ma.masked_array( data = corr, mask = mask_infobar )
        corr_infobar = np.ma.masked_array( data = corr, mask = ~ mask_infobar )
        im = ax.imshow( corr_matrix, cmap = color_map )
    if color_map_diversing : # set color_limit according to the type of a colormap
        im.set_clim( - color_limit , color_limit ) # set range of a colormap with color_limit value
    else :
        im.set_clim( 0 , color_limit ) # set range of a colormap with color_limit value
    if show_colorbar : # show a colorbar if 'show_colorbar' is True
        fig.colorbar( im, ax = ax )
    if add_Log2FC or n_factors is not None :
        if ax_infobar is None :
            im_infobar = ax.imshow( corr_infobar, cmap = 'bwr' )
        else : # when inforbar is plotted separately
            im_infobar = ax_infobar.imshow( arr_infobar, cmap = 'bwr' )
        im_infobar.set_clim( - color_limit_infobar , color_limit_infobar ) # set range of a colormap with color_limit value
        if show_colorbar : # show a colorbar if 'show_colorbar' is True
            fig.colorbar( im_infobar, ax = ax_infobar )
    if Gene_IDs_interest is not None : # if a list of Gene_IDs has been given, locate the positions of genes by overlaying a mask
        valid_Gene_IDs = list( set( df.index.values ).intersection( Gene_IDs_interest ) ) # retrive valid Gene_IDs (Gene_IDs that can be found in the dataframe)
        print( 'number (percentage) of valid Gene_IDs :', len( valid_Gene_IDs ), '(', len( valid_Gene_IDs ) / len( Gene_IDs_interest ) * 100, ')' )
        mask = deepcopy( df ) # copy dataframe to create a mask that will be used to locate a given list of genes
        mask.iloc[ :, : ] = 0
        mask[ valid_Gene_IDs ] = 1 # change colors of columns of a given list of gene_ids to white
        mask.loc[ valid_Gene_IDs ] = 1 # change colors of rows of a given list of gene_ids to white
        if add_Log2FC : # adjust the shape of a mask when Log2FC bar is added
            mask = np.vstack( ( mask.values, arr_mask ) ).T
        ax.imshow( mask, cmap = 'binary', alpha = mask_alpha ) # overlay a mask onto the correltion matrix
        title = title + '\n' + Genes_of_interest_name + ' marked by black crosses'
    elif mask_thres is not None : # if a threshold level has been given to locate the positions that has higher value than the threshold.
        mask = np.zeros( np.shape( df ) )
        mask[ df.values > mask_thres ] = 1
        if add_Log2FC or n_factors is not None :
            mask = np.vstack( ( mask.T, arr_mask ) ).T
        ax.imshow( mask, cmap = 'binary', alpha = mask_alpha ) # overlay a mask onto the correltion matrix
        title = title + '\nvalues > ' + str( round( mask_thres, 2 ) ) + ' marked by black shades'
    ax.grid( show_grid ) # hide or show grid
    if show_title :
        ax.set_title( title )
    if ax_infobar is not None :
        ax_infobar.grid( show_grid ) # hide or show grid
        ax_infobar.set_title( label_infobar )
    return ax
    #plt.tight_layout( ) # remove the grid and excessive layout to save space


# In[ ]:


"""def Calculate_correl_mat_DataFrame( df, Gene_IDs ) :
    '''    Calculate correlation matrix and Return correlation matrix dataframe    '''
    df = df[ Gene_IDs ] # reorder columns according to the clustering result
    corr = stats.spearmanr( df.values, nan_policy = 'omit' ).correlation # calculate correlation matrix DataFrame
    return pd.DataFrame( corr, index = Gene_IDs, columns = Gene_IDs ) # create DataFrame of clustered correlation matrix"""


# In[ ]:


def Find_position_in_correl_mat__Gene( df = None, Gene_of_interest = 'CRBN' ) : # set gene of interest. for example, ACAT1 as a mitochondrial marker
    '''
    Plot data of a given gene in a given clustered correlation matrix dataframe and save the plot
    '''
    gene_id = Gene_2_Gene_ID( Gene_of_interest ) # retrive gene_id
    if gene_id == -1 :
        print( 'invalid gene or gene do not exist in proteome data' )
        return
    Gene_position = np.where( df.columns.values == gene_id )[ 0 ][ 0 ] # retrive a position of the gene
    print( 'gene position :', Gene_position ) # print the position of a given gene 


# In[ ]:


def Fine_cluster__Gene_or_Position_from_tree( cluster_tree, gene = None, position = None, d = 0 ) : # set gene of interest. for example, ACAT1 as a mitochondrial marker
    '''
    Find a sub cluster of a given gene in a tree of clusters
    recursive method to visualize a dictionary-based tree
    '''
    if gene is not None :
        position = None # a given gene can overide a given position 
        gene = Gene_2_Gene_ID( gene ) # convert gene (Symbol or Gene_ID) into Gene_ID
        if gene == -1 :
            print( 'invalid gene or gene do not exist in proteome data' )
            return
    elif position is None :
        print( 'Please enter gene or position to find a specific clusters' )
        return
    start_position = cluster_tree[ 'cluster_start' ]
    end_position = cluster_tree[ 'cluster_start' ] + len( cluster_tree[ 'Gene_IDs' ] )
    if ( gene is not None and gene in cluster_tree[ 'Gene_IDs' ] ) or ( position is not None and start_position <= position and position < end_position ) :
        print( "\t" * d, cluster_tree[ 'cluster_name' ], '\tstart and end position =', start_position, ':', end_position )
        if cluster_tree[ 'subcluster_tree' ] is not None :
            for key, val in cluster_tree[ 'subcluster_tree' ].items():
                Tree__find_cluster__Gene_or_Position( val, gene = gene, position = position, d = d + 1 ) # recursion
        elif gene is not None : # if this is the last leaf node and gene is given, print the exact position of the gene 
            print( "\t" * ( d + 1 ), dict_ID_Symbol_simple[ gene ], 'position :', start_position + np.where( np.array( cluster_tree[ 'Gene_IDs' ] ) == gene )[ 0 ][ 0 ] ) # retrive a position of the gene and print the position of a given gene 


# In[ ]:


def PLOT__A_Gene_Data_Correl_Mat( df, Gene_of_interest, df_name = 'unnamed', save_fig = False, alpha = 0.5, graph_folder = None ) : # set gene of interest. for example, ACAT1 as a mitochondrial marker
    '''
    Plot data of a given gene in a given clustered correlation matrix dataframe and save the plot
    a clustered correlation matrix can be given by df argument (in this case, sample_type is used for output filename)
    '''
    Gene_ID = Gene_2_Gene_ID( Gene_of_interest )
    if Gene_ID == -1 :
        return -1 # if a given Gene_ID is invalid, print error value
    Gene_Symbol = dict_ID_Symbol_simple[ Gene_ID ]
    data = df[ int( Gene_ID ) ].values # retrive data of the given gene
    plt.plot( range( len( data ) ), data, '.', alpha = alpha ) # plot the clustered correlation scores of the given gene with the other genes 
    Gene_position = np.where( df.columns.values == Gene_ID )[ 0 ][ 0 ] # retrive a position of the gene
    plt.plot( Gene_position, 1, '.', label = Gene_Symbol ) # locate the given gene on a plot
    plt.ylabel( 'Spearman Correlation coefficient with ' + Gene_Symbol )
    plt.title( Gene_Symbol + ' entry in a clustered correlation matrix' )
    plt.legend( )
    if save_fig : # save figure if save_fig is True
        plt.savefig( graph_folder + Gene_Symbol + 'clustered_corr_matrix__' + df_name + '.png', dpi = 200 )
        plt.close( )


# In[ ]:


def EXPLORE_Cluster( df_clus, dict_tree_clus, df_data = None, df_name = '', color_map = 'bwr', add_Log2FC = False ) :
    ''' interactive method that can be used to explore clusters of 'df_clus' in 'dict_tree_clus'
    'df_data' is for drawing Log2FC T/N along with the matrix 
    it is better and faster to use %matplotlib inline instead of qt, to store all the plots that has been drawn 
    current_cluster_name can be saved and returned at the end of this method as a list '''
    current_tree_node = dict_tree_clus
    list_parent_tree_node = list( ) # a stack that will store parant nodes for trace back
    saved_cluster_names = list( )
    while True :
        current_cluster_name = current_tree_node[ 'cluster_name' ]
        current_cluster_Gene_IDs = current_tree_node[ 'Gene_IDs' ]
        print( 'n_current_cluster_Gene_IDs :', len( current_cluster_Gene_IDs ) )
        Plot_corr_matrix__Gene_IDs( df = df_clus, df_data = df_data, show_inforbar_separately = False, color_map = color_map, add_Log2FC = add_Log2FC, Gene_IDs_interest = current_cluster_Gene_IDs, Genes_of_interest_name = current_cluster_name, df_name = df_name ) # plot current cluster on relation matrix
        plt.show( )
        subclusters = current_tree_node[ 'subcluster_tree' ]
        if subclusters is None : # if subcluster_tree is None, skip displaying subcluster_tree
            print( 'current node is leaf node. No subclusters exists' )
        else :
            dict_index__clus_name__num_genes = dict( ) # build a dataframe to display subcluster_tree
            cluster_index = 0
            for subcluster_name, subcluster_node in subclusters.items( ) :
                dict_index__clus_name__num_genes[ cluster_index ] = dict( subcluster_name = subcluster_name, num_genes = len( subcluster_node[ 'Gene_IDs' ] ) )
                cluster_index += 1
            display( pd.DataFrame( dict_index__clus_name__num_genes ).transpose( ) )
        decision = '' # set an initial value of input 
        while decision == '' : # get an input
            print( "press 'q' and enter to quit exploring. press 'b' and enter to go back to parent cluster. press 's' to name of current cluster" )
            decision = input( )
            if decision == 's' : # if 's' is entered, save current name to the list
                saved_cluster_names.append( current_cluster_name ) # save current cluster name to the list, which will be returned once this interactive method ends
                print( 'current cluster_name has been saved' )
                decision = ''
        plt.close( ) # close previously opened figure
        if decision == 'q' :
            break # quit exploring clusters
        elif decision == 'b' : # go back to parent node
            if len( list_parent_tree_node ) == 0 : # if a stack list_parent_tree_node is empty, then current node is root_node, and print message and receive another input
                print( 'current cluster is the biggest cluster and do not have parent cluster')
                continue
            else :
                current_tree_node = list_parent_tree_node.pop( ) # go back to immediate parent_tree_node
                continue 
        elif ( subclusters is not None ) and int( decision ) in dict_index__clus_name__num_genes : # if subcluster_tree is not empty, change current cluster to a subcluster according to the given index
            next_cluster_name = dict_index__clus_name__num_genes[ int( decision ) ][ 'subcluster_name' ] # retrive name of next cluster
            list_parent_tree_node.append( current_tree_node ) # push current node to the the parent node stack
            current_tree_node = subclusters[ next_cluster_name ] # set current cluster node
    print( 'a list of saved cluster names returned')
    return saved_cluster_names # return a list of saved cluster names


# In[ ]:


def PLOT_magnigfied_correl_matrix( df, df_data, start, end, show_grid = False, color_limit = 1, size = 'auto' ) :
    ''' For clustered relation matrix (let's call ratio matrix or correlation matrix all together as relation matrix), 
        Magnify the matrix from 'start' to 'end', and Add Log2FC data at the rightmost column (df_data should be given through an argument)
        Also, return a DataFrame with indices, Gene_Name_Symbols, and Log2FC so that Two-gene plot can be subsequently drawn '''
    end += 1
    n_genes = end - start # retrive number of genes
    Gene_IDs = df.index.values[ start : end ] # retrive Gene_IDs
    Gene_Name_Symbols = List_Gene_ID__2__List_Gene_Symbol( Gene_IDs, add_gene_name = True ) # retrive Gene_Name(Gene_Symbol) annotations
    if size == 'auto' : # automatically set size of the plot
        size = n_genes / 3
    fig, ax = plt.subplots( figsize = ( size, size ) )
    Log2FC = Calculate_Log2FC_p_value__A_vs_B( df_data, Gene_IDs = Gene_IDs ).loc[ Gene_IDs, 'Log2_Fold_Change' ].values # calculate Log2FC
    cax = ax.matshow( np.vstack( ( df.loc[ Gene_IDs, Gene_IDs ].values.T, Log2FC ) ).T, cmap = 'bwr', vmin = - color_limit, vmax = color_limit )
    plt.xticks( range ( len( Gene_Name_Symbols ) + 1 ), Gene_Name_Symbols + [ 'Log2FC T/N' ], rotation = 270 ); # add column labels
    plt.yticks( range ( len( Gene_Name_Symbols ) ), Gene_Name_Symbols );
    # Add the colorbar legend
    cbar = fig.colorbar( cax, shrink = 0.8 )
    plt.grid( show_grid )
    return pd.DataFrame( dict( Gene_Name_Symbols = Gene_Name_Symbols, Log2FC_T_vs_N = Log2FC ), index = np.arange( start, end ) )


# In[ ]:


def CLUSTER_GET_Gene_IDs_by_cluster_name( dict_tree_clus, cluster_name ) :
    ''' Recursive method to retrive Gene_IDs of a given name of cluster. return -1 if no such cluster_name exist '''
    if dict_tree_clus[ 'cluster_name' ] == cluster_name : # if current cluster_name match a given cluster_name, return Gene_IDs of it
        return dict_tree_clus[ 'Gene_IDs' ]
    if dict_tree_clus[ 'subcluster_tree' ] is None : # if current cluster_name is not matched with a given name and a leaf node, return -1
        return - 1 
    else :
        for subcluster_node in dict_tree_clus[ 'subcluster_tree' ].values( ) : # for each subcluster, recursively search the cluster name 
            search_result = CLUSTER_GET_Gene_IDs_by_cluster_name( subcluster_node, cluster_name = cluster_name ) # recursion
            if search_result != -1 :
                return search_result # if a matched cluster has been found, return the results, assuming that there is only one unique name
    return - 1


# ### Utility functions for Matplotlib Plotting

# In[ ]:


def MATPLOTLIB_savefig( title, dpi = 200, folder = None, close_fig = True, format = '.png' ) :
    if '.' not in format :
        format = '.' + format
    plt.savefig( folder + To_window_path_compatible_str( title ) + format, dpi = 200, bbox_inches = 'tight' ) # prevent x or y labels from cutting off
    if close_fig :
        plt.close( )


# In[ ]:

def MPL_SAVE_svg_png_pdf( fig_name, ** dict_save_fig ) :
    ''' With the given 'fig_name', save fiqures in both svg and png format '''
    MATPLOTLIB_savefig( fig_name, format = '.svg', close_fig = False, ** dict_save_fig )
    MATPLOTLIB_savefig( fig_name, format = '.pdf', close_fig = False, ** dict_save_fig )
    MATPLOTLIB_savefig( fig_name, ** dict_save_fig )
MPL_SAVE_svg_png = MPL_SAVE_svg_png_pdf

# In[ ]:


def MATPLOTLIB_basic_configuration( font_size = None, font_size_axes_title = None, font_size_axes_label = None, font_size_xtick = None, font_size_ytick = None, font_size_legend = None, font_size_figure_title = None, x_label = None, y_label = None, title = None, x_scale = None, y_scale = None, show_grid = True, show_legend = False, savefig = False, y_lim = None, x_lim = None, save_file_name = None, folder = None, format = '.png', show_colorbar = False ) :
    ''' A basic function for confiquring a matplotlib plot '''
    # set font sizes
    if font_size is not None : plt.rc( 'font', size = 20 ) # controls default text sizes
    if font_size_axes_title is not None : plt.rc( 'axes', titlesize = 20 ) # fontsize of the axes title   
    if font_size_axes_label is not None : plt.rc( 'axes', labelsize = 20 ) # fontsize of the x and y labels
    if font_size_xtick is not None : plt.rc( 'xtick', labelsize = 20 ) # fontsize of the x tick labels  
    if font_size_ytick is not None : plt.rc( 'ytick', labelsize = 20 ) # fontsize of the y tick labels  
    if font_size_legend is not None : plt.rc( 'legend', fontsize = 20 ) # legend fontsize              
    if font_size_figure_title is not None : plt.rc( 'figure', titlesize = 50 ) # fontsize of the figure title  
    if x_label is not None : plt.xlabel( x_label )
    if y_label is not None : plt.ylabel( y_label )
    if title is not None : plt.title( title )
    if x_scale is not None : plt.xscale( x_scale )
    if y_scale is not None : plt.yscale( y_scale )
    if x_lim is not None :
        if isinstance( x_lim, ( tuple, list ) ) : plt.xlim( left = x_lim[ 0 ], right = x_lim[ 1 ] )
        elif isinstance( x_lim, dict ) : plt.xlim( **x_lim )
    if y_lim is not None :
        if isinstance( y_lim, ( tuple, list ) ) : plt.ylim( bottom = y_lim[ 0 ], top = y_lim[ 1 ] )
        elif isinstance( y_lim, dict ) : plt.ylim( **y_lim )
    plt.grid( show_grid )
    if show_legend : plt.legend( )
    if savefig :
        if save_file_name is None : # if 'save_file_name' is not given 
            if title is None : title = 'Unnamed Plot_' + TIME_GET_timestamp( ) # if title is not given, put a default title to save a plot
            MATPLOTLIB_savefig( title = title, folder = folder, format = format )
        else : MATPLOTLIB_savefig( title = save_file_name, folder = folder, format = format )
    if show_colorbar : plt.colorbar( )
MPL_basic_configuration = MATPLOTLIB_basic_configuration


# ### Functions for drawing basic plots using Matplotlib 

# ##### Functions for matplotlib setting 

# In[ ]:


def SETTING_MPL_DPI_for_INLINE_plot( dpi = None ) :
    '''   Set resolution of an inline plot   '''
    dpi = 120 if dpi is None else dpi # set default dpi
    mpl.rcParams[ 'figure.dpi' ] = dpi
MPLSETTING_dpi_for_INLINE_plot = SETTING_MPL_DPI_for_INLINE_plot


# In[ ]:


l_setting_mpl_default = cycler( color = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ] )
def SETTING_MPL_Default_Property_Cycle( ** dict_property ) :
    '''   Set default property cycle for matplotlib. default color property cycle is 'bgrcmyk'   '''
    l_default_color = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ]
    mpl.rcParams[ 'axes.prop_cycle' ] = cycler( color = l_default_color, ** dict_property ) if dict_property.get( 'color', None ) is None else cycler( ** dict_property )


# In[1]:


def MPL_1D_Sort_Plot( data, figsize = ( 5, 3.5 ), annotate_xy_pos_first_column_label = ( 0.05, 1.09 ), color_alpha = 0.5, color_threshold = 0, line_stype = '.-', x_label = 'Sorted Entries', title = '', savefig = False, color_above = 'g', color_below = 'r', color_percentile_alpha = 0.5, color_percentile_thres = None, color_percentile_lower = 'b', color_percentile_upper = 'orange', thres_n_points = 10000, ** dict_mpl_basic_configure ) :
    ''' (1) Convert iterable data like series or list into np.ndarray using 'TYPE_Convert_NP_Array' (2) Sort, (3) Visualize on a plot using green and red colors 
    to visualize deviation from the given threshold. if color_percentile_thres is not None, annotate upper and lower percentiles with the color given by color_percentile arguments
    if 'data' is pandas.DataFrame with two columns, first sort values in the second column by the first column, visualize, and annotate unique_entries of the first column on the plot. The NaN values in the first column will be ignored. '''
    bool_flag_sort_using_two_columns = isinstance( data, pd.DataFrame ) and len( data.columns.values ) == 2 # 'sort_using_two_columns' if a DataFrame with two columns are given as a 'data'.
    if isinstance( data, ( pd.DataFrame, pd.Series ) ) : # set default title and y_label by using pandas.Series name if not given
        data_name = data.name if isinstance( data, pd.Series ) else data.columns.values[ 1 ]
        if data_name : # if data_name is not None
            if 'y_label' not in dict_mpl_basic_configure : dict_mpl_basic_configure[ 'y_label' ] = data.name if isinstance( data, pd.Series ) else data.columns.values[ 1 ]
            if len( title ) == 0 : title = data_name # set data_name as default title for the plot
    if not bool_flag_sort_using_two_columns : # convert data as an numpy array and sort the data if 'series__sort_by_index_first' is set to False or given 'data' is not pandas.Series
        data = TYPE_Convert_NP_Array( data, dtype = float ) 
        if type( data ) is not np.ndarray : return -1
    if bool_flag_sort_using_two_columns : data = data.dropna( ) # remove np.nan values if present
    else :
        arr_mask_isnan = np.isnan( data )
        if arr_mask_isnan.any( ) : data = data[ ~ arr_mask_isnan ] 
    int_ratio_shrinkage = int( len( data ) / thres_n_points ) if len( data ) > thres_n_points else 1 # shrink the number of data values before sorting for efficient plotting. # limit the number of points beling plotting below 'thres_n_points'
    if int_ratio_shrinkage > 1 : 
        print( "Since len( data ) = {} + ({} NaN values) > thres_n_points = {}, Perform random sampling of the data (one value for every {} values) for efficient visualization".format( len( data ), arr_mask_isnan.sum( ), thres_n_points, int_ratio_shrinkage ) )
        data = data.iloc[ : : int_ratio_shrinkage ] if bool_flag_sort_using_two_columns else data[ : : int_ratio_shrinkage ] 
    data = data.sort_values( list( data.columns ), inplace = False, ignore_index = True ) if bool_flag_sort_using_two_columns else np.sort( data ) # sort data
    arr_data = data.iloc[ :, 1 ].values if bool_flag_sort_using_two_columns else data # plot data
    x_range = np.arange( len( arr_data ) )
    x_axis = np.full_like( x_range, color_threshold )
    fig, ax = plt.subplots( 1, 1, figsize = figsize )
    ax.plot( arr_data, line_stype, color = 'k' ) 
    if color_percentile_thres is not None : # fill colors in the plot
        index_thres = int( float( len( arr_data ) ) * color_percentile_thres / 100 )
        x_axis_plus_max = np.full_like( x_range, arr_data[ - 1 ] )
        ax.fill_between( x_range[ : index_thres ], x_axis_plus_max[ : index_thres ], x_axis[ : index_thres ], facecolor = color_percentile_lower, interpolate = True, alpha = color_percentile_alpha )
        ax.fill_between( x_range[ - index_thres : ], x_axis_plus_max[ - index_thres : ], x_axis[ - index_thres : ], facecolor = color_percentile_upper, interpolate = True, alpha = color_percentile_alpha )
    ax.fill_between( x_range, arr_data, x_axis, where = arr_data >= x_axis, facecolor = color_above, interpolate = True, alpha = color_alpha )
    ax.fill_between( x_range, arr_data, x_axis, where = arr_data <= x_axis, facecolor = color_below, interpolate = True, alpha = color_alpha )
    plt.sca( ax ) # set x_ticks properly after shrinkage
    arr_xticks = plt.xticks( )[ 0 ][ 1 : -1 ]
    plt.xticks( arr_xticks, ( arr_xticks * int_ratio_shrinkage ).astype( int ) )
    if bool_flag_sort_using_two_columns : # annotate unique entries in the first columns by which values of second columns were first sorted.
        l_unqiue_entry = sorted( data.iloc[ :, 0 ].unique( ) )
        dict_unique_entry_to_int_representation = dict( ( unique_entry, int_representation ) for int_representation, unique_entry in enumerate( l_unqiue_entry ) )
        dict_int_representation_to_unique_entry = dict( ( int_representation, unique_entry ) for int_representation, unique_entry in enumerate( l_unqiue_entry ) )
        data.iloc[ :, 0 ] = list( dict_unique_entry_to_int_representation[ entry ] for entry in data.iloc[ :, 0 ].values )
        l_start_of_unique_entry = [ 0 ] + list( np.where( np.diff( data.iloc[ :, 0 ].values ) )[ 0 ] ) + [ len( data ) ]
        for int_representation, unique_entry in enumerate( l_unqiue_entry ) :
            x_pos = ( l_start_of_unique_entry[ int_representation ] + l_start_of_unique_entry[ int_representation + 1 ] ) / 2
            ax.annotate( unique_entry, xy = ( x_pos, 1.02 ), xycoords = ( "data", "axes fraction" ), ha = "center" )
        ax.annotate( data.columns.values[ 0 ], xy = annotate_xy_pos_first_column_label, xycoords = ( "axes fraction", "axes fraction" ), ha = "center" )
    MATPLOTLIB_basic_configuration( x_label = x_label, title = TIME_GET_timestamp( ) + '\n' + title, savefig = savefig, ** dict_mpl_basic_configure )


# In[ ]:


def MPL_Internal_Util_get_min_max_for_plotting_an_axis( arr ) :
    min_value, max_value = arr.min( ), arr.max( )
    if min_value * max_value > 0 :
        if min_value > 0 :
            min_value = 0
        else :
            max_value = 0
    return min_value, max_value


# In[ ]:


def MPL_Scatter_Align_Two_Series( s_1, s_2, ls = '', marker = 'o', alpha = 0.5, show_xy_axis = ( True, True ), annotate_xy_axis = ( False, False ), annotate_color = 'blue', annotate_labels = None, annotate_outlier_in_both_x_and_y = True, annotate_distance_factor = 1, n_std_for_outliers = 3, label_split_char = '_', figsize = ( 5, 5 ), ** dict_mpl_basic_config ) :
    ''' Align two series and plot a scatter plots. Annotate outliers by default. Annotations of list of points are also availble via 'annotate_labels' arguments. '''
    s_1, s_2 = s_1.dropna( ).align( s_2.dropna( ), join = 'inner' ) # align two series
    arr_1, arr_2, arr_labels = s_1.values.astype( float ), s_2.values.astype( float ), s_1.index.values
    if label_split_char is not None and type( arr_labels[ 0 ] ) is str :
        arr_labels = np.array( list( label.split( '_' )[ 0 ] for label in arr_labels ) ) # split labels to reduce complexity of annotations
    if not isinstance( arr_labels[ 0 ], ( str, np.str_ ) ) :
        labels, mask = List_Gene_ID__2__List_Gene_Symbol( arr_labels, return_mask_mapped = True )
        arr_labels = arr_labels.astype( object )
        arr_labels[ mask ] = labels
        arr_labels[ ~ mask ] = 'Gene Not Mapped'
    min_1, max_1 = MPL_Internal_Util_get_min_max_for_plotting_an_axis( arr_1 ) # retrive min, max, range of x and y values for plotting
    min_2, max_2 = MPL_Internal_Util_get_min_max_for_plotting_an_axis( arr_2 )
    range_x, range_y = max_1 - min_1, max_2 - min_2
    fig, ax = plt.subplots( 1, 1, figsize = figsize )
    if show_xy_axis is not None : # plot x and y axis
        step_x, step_y = range_x * 0.05, range_y * 0.05 # calculate steps for x and y axes
        arr_x_range, arr_y_range = np.arange( min_1 - step_x, max_1 + step_x + range_x * 0.01, step_x ), np.arange( min_2 - step_y, max_2 + step_y + range_y * 0.01, step_y )
        if show_xy_axis[ 0 ] :
            ax.plot( arr_x_range, np.zeros_like( arr_x_range ), color = 'black', lw = 1 )
        if show_xy_axis[ 1 ] :
            ax.plot( np.zeros_like( arr_y_range ), arr_y_range, color = 'black', lw = 1 )
    ax.plot( arr_1, arr_2, ls = ls, alpha = alpha, marker = marker ) # plot all data points in the two aligned series
    if annotate_xy_axis is not None and annotate_xy_axis[ 0 ] or annotate_xy_axis[ 1 ] : # annotate points if valid 'annotate_xy_axis' is given
        if annotate_labels is None : # if labels for annotations are not given, retrive annotate outliers by default 
            mask_arr_1, mask_arr_2 = OUTLIERS_GET_mask_for_outliers( np.vstack( ( arr_1, arr_2 ) ), n_std_for_outliers = n_std_for_outliers ) # retrive mask for outliers
            if annotate_xy_axis[ 0 ] and annotate_xy_axis[ 1 ] : # retrive a mask for annotation of outliers for both x and y axis 
                mask_arr = mask_arr_1 & mask_arr_2 if annotate_outlier_in_both_x_and_y else mask_arr_1 | mask_arr_2
            else : # retrive a mask for annotation of outliers for either x or y axis 
                mask_arr = mask_arr_1 if annotate_xy_axis[ 0 ] else mask_arr_2
        else :
            mask_arr = GET_MASK_of_intersection( arr_labels, annotate_labels )
        ax.plot( arr_1[ mask_arr ], arr_2[ mask_arr ], 'o', color = annotate_color, alpha = 0.7 ) # annotate outliers
        for x, y, label in zip( arr_1[ mask_arr ], arr_2[ mask_arr ], arr_labels[ mask_arr ] ) :
            ax.annotate( label, ( x, y ), xytext = ( x + range_x * 0.01 * annotate_distance_factor, y + range_y * 0.01 * annotate_distance_factor ) )
    if 'x_label' not in dict_mpl_basic_config : dict_mpl_basic_config[ 'x_label' ] = s_1.name # use pandas Series names as default axis labels
    if 'y_label' not in dict_mpl_basic_config : dict_mpl_basic_config[ 'y_label' ] = s_2.name
    MATPLOTLIB_basic_configuration( ** dict_mpl_basic_config )


# In[ ]:


def MPL_2D_Hist( s_1, s_2, colorbar = True, logscale = True, bins = 100, figsize = ( 7, 5 ), cmap = 'inferno', ** dict_mpl_basic_config ) : # 2020-07-23 23:34:28 
    ''' Align two series and plot a scatter plots. Annotate outliers by default. Annotations of list of points are also availble via 'annotate_labels' arguments. '''
    if isinstance( s_1, ( pd.Series ) ) and isinstance( s_2, ( pd.Series ) ) : # align two series if series is given
        s_1, s_2 = s_1.dropna( ).align( s_2.dropna( ), join = 'inner' ) 
        arr_1, arr_2 = s_1.values.astype( float ), s_2.values.astype( float )
        if 'x_label' not in dict_mpl_basic_config and s_1.name is not None : dict_mpl_basic_config[ 'x_label' ] = s_1.name # use series name as default axis label
        if 'y_label' not in dict_mpl_basic_config and s_2.name is not None : dict_mpl_basic_config[ 'y_label' ] = s_2.name
    else : arr_1, arr_2 = s_1, s_2 # if numpy arrays were given
    fig, ax = plt.subplots( 1, 1, figsize = figsize )
    colors_norm = colors.LogNorm( ) if logscale else colors.Normalize( ) # set normalization method based on 'logscale' argument
    h, xedges, yedges, im = ax.hist2d( arr_1, arr_2, bins = bins, norm = colors_norm, cmap = cmap )
    if colorbar : plt.colorbar( im, ax = ax ) # show colorbar if 'colorbar' is true
    MATPLOTLIB_basic_configuration( ** dict_mpl_basic_config )


# In[ ]:


def PLOT_Scatter_Annotation( arr_coordinates, arr_labels, cmap = None, alpha = 0.7, figsize = ( 8.0, 6.0 ), marker = 'o', dict_ax_scatter = dict( ), ** dict_matplotlib_basic_setting ) :
    '''  'arr_coordinates' should be numpy array with 2 columns, x and y coordinates, and 'arr_labels' should be a list-like object containing labels. 
    If arr_labels contains np.nan and it is pd.Series, remove np.nan labels from the given arr_labels and entries with np.nan labels from arr_coordinates  '''
    if NUMPY_UTIL_Does_contain_NaN( arr_labels ) : # if arr_labels contains np.nan and it is pd.Series, remove np.nan labels from the given arr_labels and entries with np.nan labels from arr_coordinates
        arr_coordinates = arr_coordinates[ ~ pd.isnull( arr_labels ) ]
        arr_labels = arr_labels.dropna( )
    arr_labels = arr_labels if isinstance( arr_labels, ( np.ndarray ) ) else np.array( arr_labels, dtype = object ) # convert arr_labels to a numpy object
    s_label = LIST_COUNT( arr_labels, duplicate_filter = None ).sort_values( ascending = False )
    fig, ax = plt.subplots( figsize = figsize )  
    l_colors = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ] if cmap is None else COLORMAP_GET_list_of_colors_from_list_of_values( np.arange( len( s_label ) ), return_RGB_hex_str = True, cmap = cmap ) # set color map for plotting 
    for label, count, color in zip( s_label.index, s_label.values, l_colors ) : # draw markers for samples in each cluster
        x_coordinates, y_coordinates = arr_coordinates[ arr_labels == label ].T
        ax.scatter( x_coordinates, y_coordinates, label = '{Label} (N={N_labels})'.format( Label = label, N_labels = int( count ) ), alpha = alpha, color = color, marker = marker, ** dict_ax_scatter ) # for each cluster, plot markers with different colors
    plt.legend( )
    MATPLOTLIB_basic_configuration( ** dict_matplotlib_basic_setting )


# ##### functions utilizing imshow

# In[ ]:


def IMSHOW_expand_and_visualiz_narrow_arr( df, label_columns = None, width_column = None, cmap = 'bwr', show_grid = False, color_limits = ( -1, 1 ), gene_list_annotation = None, font_size_col_index = 10, font_size_col = 10 ) :
    '''  'gene_list_annotation' : See 'GENE_Annotate_df_for_plot'  '''
    if label_columns is None : # set default column labels
        label_columns = df.columns.values
    if width_column == 1 :
        arr_expanded = df.values
    else :
        arr_expanded = NUMPY_expand_narrow_arr( df, width_column = width_column ) # expand df to help visualization of DataFrame
    mask_nan = np.isnan( arr_expanded ) # retrive mask of NaN values
    width_column = int( arr_expanded.shape[ 1 ] / df.values.shape[ 1 ] ) # retrive the value of width_columns used to expand array
    fig, ax = plt.subplots( 1, 1 )
    im = ax.imshow( np.ma.masked_array( arr_expanded, mask_nan ), cmap = cmap ) # show heatmap of data
    im.set_clim( color_limits[ 0 ], color_limits[ 1 ] )
    im_nan = ax.imshow( np.ma.masked_array( np.zeros_like( arr_expanded ), ~ mask_nan ), cmap = 'Greys_r' ) # show invalid values using a black color
    im_nan.set_clim( 0, 1 )
    ax.grid( show_grid )
    Symbols_annotation, list_positions_annotation = GENE_Annotate_df_for_plot( df, gene_list_annotation ) # retrive position and symbols of genes (indices) that will be annotated
    if Symbols_annotation != -1 : # if annotation exists
        ax.set_yticks( list_positions_annotation ) # annotate index and columns 
        ax.set_yticklabels( Symbols_annotation, fontsize = font_size_col_index )
    ax.set_xticks( np.arange( int( width_column / 2 ), arr_expanded.shape[ 1 ], width_column ) )
    ax.set_xticklabels( label_columns, fontsize = font_size_col )
    ax.tick_params( top = True, bottom = False, labeltop = True, labelbottom = False ) # Let the horizontal axes labeling appear on top.
    plt.setp( ax.get_xticklabels( ), rotation = 30, ha = "left", rotation_mode = "anchor" ) # Rotate the tick labels and set their alignment.


# ### Utility Functions for Bokeh interactive Plotting

# In[ ]:


def OUTLIERS_GET_data_values_without_outliers( arr_values, n_std_for_outliers = 2 ) :
    ''' remove outliers from a 1-d array 'arr_values' (outside std * 'n_std_for_outliers') and return data_values after outlier removal '''
    mask_outliers = OUTLIERS_GET_mask_for_outliers( np.array( [ arr_values ], dtype = float ), n_std_for_outliers = n_std_for_outliers )[ 0 ] # remove outliers from list_values to map scalar values to color effectively 
    return arr_values[ ~ mask_outliers ] # retrive data values without outliers


# In[ ]:


def OUTLIERS_GET_min_max_after_outlier_removal( arr_values, n_std_for_outliers = 2 ) :
    ''' remove outliers from a 1-d array 'arr_values' (outside std * 'n_std_for_outliers') and return min and max after outlier removal '''
    if type( arr_values ) is not np.ndarray :
        arr_values = np.array( arr_values, dtype = float )
    arr_values_without_outliers = OUTLIERS_GET_data_values_without_outliers( arr_values, n_std_for_outliers = n_std_for_outliers ) # remove outliers from list_values to map scalar values to color effectively 
    return arr_values_without_outliers.min( ), arr_values_without_outliers.max( ) # retrive min and max of values without outliers 


# In[ ]:


def BOKEH_transform_arr_values_into_a_range_from_0_to_1( arr_values, n_std_for_outliers = 2, log_transform = False, inverse = False ) :
    ''' remove outliers from a 1-d array 'arr_values', retrive max and min of values without outliers, and map values to color using a given colormap and return np.array of colors. if list_values is p_values, set 'log_transform' to True to get appropriate coverage of the color '''
    if type( arr_values ) is not np.ndarray :
        arr_values = np.array( arr_values, dtype = float )
    if log_transform : # if 'log_transform' is True, transform the data into log-values. 
        arr_values = np.log( arr_values )
    if inverse : # if inverse is True, set the maximum value to 0 and change the sign so that ranks of data_values becomes the opposite while all data_values remain positive
        arr_values = arr_values.max( ) - arr_values
    else :
        arr_values = arr_values - arr_values.min( ) # set the minimun value to 0 so that all data_values become positive
    if n_std_for_outliers is None : # if None is given through 'n_std_for_outliers' argument, do not exclude outliers when retriving min and max values  
        arr_values_without_outliers_min, arr_values_without_outliers_max = arr_values.min( ), arr_values.max( )
    else :
        arr_values_without_outliers_min, arr_values_without_outliers_max = OUTLIERS_GET_min_max_after_outlier_removal( arr_values, n_std_for_outliers = n_std_for_outliers ) # retrive min and max of values without outliers 
    arr_values_transformed = ( arr_values - arr_values_without_outliers_min ) / arr_values_without_outliers_max
    arr_values_transformed[ arr_values_transformed > 1 ] = 1 # set transformed values of outliers to maxinum (1) and minimum (0) values 
    arr_values_transformed[ arr_values_transformed < 0 ] = 0
    return arr_values_transformed


# In[ ]:


def COLORMAP_GET_list_of_colors_from_list_of_values( list_values, cmap = 'viridis', n_std_for_outliers = 2, return_RGB_hex_str = False, log_transform = False, diverging_colormap = False ) :
    ''' remove outliers fron list_values, retrive max and min of values without outliers, and map values to color using a given colormap and return np.array of colors. if list_values is p_values, set 'log_transform' to True to get appropriate coverage of the color '''
    if type( list_values ) is not np.ndarray :
        list_values = np.array( list_values, dtype = float )
    if log_transform : # if 'log_transform' is True, transform the data into log-values. 
        list_values = np.log( list_values )
    if n_std_for_outliers is None : # if 'n_std_for_outliers' is None, do not exclude outliers when retriving min and max values
        arr_values_without_outliers_min, arr_values_without_outliers_max = list_values.min( ), list_values.max( )
    else :
        arr_values_without_outliers_min, arr_values_without_outliers_max = OUTLIERS_GET_min_max_after_outlier_removal( list_values, n_std_for_outliers = n_std_for_outliers ) # retrive min and max of values without outliers 
    if arr_values_without_outliers_max == 0 : # if resulting 'arr_values_without_outliers_max' is zero, retrive maximum value without considering zero values
        arr_values_without_outliers_min, arr_values_without_outliers_max = OUTLIERS_GET_min_max_after_outlier_removal( list_values[ list_values != 0 ], n_std_for_outliers = n_std_for_outliers )
    if cmap in cm.__dict__ : # if given color_map is valid, use the scalar color mapper of a given color_map name
        colormapper = cm.__dict__[ cmap ]
    else : # if the given color_map is invalid, return an error message and end the method
        print( 'invalid color_map' )
        return -1 
    if diverging_colormap :
        limit_val = np.abs( [ arr_values_without_outliers_min, arr_values_without_outliers_max ] ).max( )
        arr_values_without_outliers_min, arr_values_without_outliers_max = - limit_val, limit_val
    if arr_values_without_outliers_min < 0 : # if minimum value is negative (regardless of the sign of arr_values_without_outliers_max), increase arr_values_without_outliers_max by arr_values_without_outliers_min so that arr_values_without_outliers_max becomes 1 after transformation
        arr_values_without_outliers_max -= arr_values_without_outliers_min
    arr_color = np.array( list( colormapper( int( 255 * float( value - arr_values_without_outliers_min ) / arr_values_without_outliers_max ) )[ : 3 ] for value in list_values ) ) # map a given list of values with color_map of a given color_map name to colors and return the list of colors as a np.array
    if return_RGB_hex_str : # if 'return_RGB_hex_str' is set to True, convert RGB color tuple into RGB hex string and return list of RGB hex strings
        arr_color_0_to_255 = np.array( arr_color * 255, dtype = int ) # convert 0 ~ 1 RGB values into 0 ~ 255 RGB integer values
        return np.array( list( '#%02x%02x%02x' % ( color_0_to_255[ 0 ], color_0_to_255[ 1 ], color_0_to_255[ 2 ] ) for color_0_to_255 in arr_color_0_to_255 ) ) # convert RGB color tuple into RGB hex string
    else : # if 'return_RGB_hex_str' is set to False, return list of RGB color tuples
        return arr_color


# ### Plotting Functions for Exploratory Analysis (Normalization)

# In[ ]:


def VERIFICATION_boxplot_distribution_of_ranks_of_genes_for_each_sample( df, sample_list, Gene_Set = None, save_fig = False, show_plot = True, name_df =  'df_proteome_unshared', name_Gene_Set = 'Unnamed', name_sample_list = 'Unnamed' ) :
    ''' Calculate Rank of data of samples for each gene, and calculate average rank for each sample, and draw a distribution of ranks of genes in a given set of genes for each sample as a boxplot 
    Return series containing data of arr_sample_ranks_median with sample_list as indices  '''
    if Gene_Set is not None :
        df = PANDAS_Subset( df, Gene_Set )
    df = df[ sample_list ]
    arr_sample_ranks = df.values.argsort( axis = 1 ).argsort( axis = 1 ) / len( sample_list ) * 100 # normalize Rank into range from 0 to 100 to aid visualization
    arr_sample_ranks_median = np.median( arr_sample_ranks, axis = 0 )
    indices_sorting_rank_median = arr_sample_ranks_median.argsort( )
    sample_list_sorted, _ = ARRAY_sort_one_list_based_on_another( sample_list, arr_sample_ranks_median )
    if show_plot : # if show_plot is True, draw a plot of boxplots
        plt.figure( figsize = ( 18, 12 ) ) 
        plt.boxplot( arr_sample_ranks.T[ indices_sorting_rank_median ].T, widths = 0.85, showmeans = True, notch = True, showfliers = False, whis = 0.25 ) # draw box plot
        title = dataset_name + " Distribution of Ranks of '{}' Genes in {} Among {} Samples".format( name_Gene_Set, name_df, name_sample_list )
        plt.title( title )
        plt.ylabel( 'Distribution of Ranks of Given Genes for Each Sample.' )
        plt.xticks( np.arange( 1, len( sample_list ) + 1 ), sample_list_sorted, rotation = 90, fontsize = 9 )
        plt.ylim( 0, 100 )
        if save_fig :
            MATPLOTLIB_savefig( title )
    return pd.Series( arr_sample_ranks_median, index = sample_list )


# ### Individual sample, Gene-level, Not interactive Basic Plotting 

# In[ ]:


def PLOT_GEL_BLOTs( sample_list, Gene_Set, df = None, method_data_heatmap = 'max', sample_labels = None, add_Gene_Name = True, width_blot = 0.8, cmap = 'Greys_r', textcolors = [ "black", "white" ], text_color_threshold = 50, threshold_bothends = None, title = None, show_colorbar = False, valfmt = "{x:.2e}", aspect = 0.3, grid_color = 'k', grid_width = 0.7, data_label = 'intensity', fontsize_col = 12, fontsize_index = 12, fontsize_blot = 12 ) :
    ''' Plot Grey_scale heatmap of data of a given DataFrame with a given sample_list and Gene_Set, as if showing a image of gel blots (for example, Western Blot) 'data_label' is a label for colorbar label
    'method_data_heatmap' : method used to set colors of heatmap.  'max' = max value become the highest value and 0 becomes lowest value. 'std_positive' = log-average + 2 * std becomes max color
            and 0 becomes lowest value. 'std_zero' = 2 * std and -2 *std becomes max and min color. A tuple of color limits can be alternatively given when data_heatmap == data '''
    if df is None :
        df = df_rna
    valid_sample_list = np.array( sample_list )[ GET_MASK_of_intersection( sample_list, df.columns.values ) ] # retrive valid sample_list while preserving the order of a given list of samples
    df = df[ valid_sample_list ] # retrive data of valid_sample_list
    df = df.replace( 0, np.nan ).dropna( ) # if mask_nan_color is None else df.replace( 0, np.nan ) # drop entries with zero values if None values are not masked. Mask zero values with np.nan values
    if isinstance( df.index.values[ 0 ], ( int, float, np.float64, np.int64 ) ) : # if the given df's indices are gene_ids
        Gene_IDs = List_Gene__2__List_Gene_ID( Gene_Set )
        valid_Gene_IDs = Gene_IDs[ GET_MASK_of_intersection( Gene_IDs, df.index.values ) ] # retrive valid Gene_IDs while preserving the order of a given list of genes
        Gene_Symbols = List_Gene_ID__2__List_Gene_Symbol( valid_Gene_IDs, add_Gene_Name )
    else : # if the given df's indices are not gene_ids
        valid_Gene_IDs, Gene_Symbols, add_Gene_Name, width_blot, fontsize_index = Gene_Set, Gene_Set, False, width_blot / 2, 17 # set values for plotting df whose indices are not gene_id 
    data = df.loc[ valid_Gene_IDs ].values
    if np.sum( data < 0 ) > 0 : # if data contain negative values, set color limit automatically
        method_data_heatmap, valfmt = ( -1 , 0 ), "{x:.2f}"
    if method_data_heatmap == 'max' :
        data_heatmap, color_limits, text_color_threshold = ( data.T / data.max( axis = 1 ) ).T * 100, ( 0, 100 ), 50 # Normalize data with its maximum values, set color limits, and threshold for text color. 
    elif method_data_heatmap == 'std_positive' :
        data_heatmap, color_limits, text_color_threshold = ( data.T / np.power( 2, np.log2( data ).mean( axis = 1 ) ) + data.std( axis = 1 ) * 2 ).T * 100, ( 0, 100 ), 50 # Normalize data with its log-avg + 2 std, set color limits, and threshold for text color. 
    elif method_data_heatmap == 'std_zero' :
        data_heatmap, color_limits, text_color_threshold = ( data.T / data.std( axis = 1 ) * 2 ).T * 100, ( -100, 100 ), 0 # Normalize data with 2 * std, set color limits, and threshold for text color. 
    elif type( method_data_heatmap ) is tuple : 
        data_heatmap, color_limits, text_color_threshold = data, method_data_heatmap, np.mean( method_data_heatmap ) # do not normalize data, set color limits, and threshold for text color. 
    else :
        return - 1
    fig, ax = plt.subplots( figsize = ( 5 + int( len( valid_sample_list ) * width_blot ) + 2 * int( add_Gene_Name ), 2 + int( len( valid_Gene_IDs ) / 2 ) + len( valid_sample_list[ 0 ] ) / 5 ) )
    ax.grid( False ) # Turn off grid
    if sample_labels is None : # draw heatmap # if sample_labels has not be separately given, use sample names to annotate each sample
        im, cbar = heatmap( data_heatmap, row_labels = Gene_Symbols, col_labels = valid_sample_list, ax = ax, show_colorbar = show_colorbar, cmap = cmap, cbarlabel = data_label, color_limit = color_limits, grid_color = grid_color, grid_width = grid_width, cbar_kw = dict( aspect = 3 ), fontsize_x = fontsize_col, fontsize_y = fontsize_index )
    else : 
        im, cbar = heatmap( data_heatmap, row_labels = Gene_Symbols, col_labels = sample_labels, ax = ax, show_colorbar = show_colorbar, cmap = cmap, cbarlabel = data_label, color_limit = color_limits, grid_color = grid_color, grid_width = grid_width, cbar_kw = dict( aspect = 3 ), fontsize_x = fontsize_col, fontsize_y = fontsize_index )
    
    if threshold_bothends is not None :
        texts = annotate_heatmap(im, data, valfmt = valfmt, textcolors = textcolors[ : : -1 ], threshold_bothends = threshold_bothends, fontsize = fontsize_blot )
    else :
        texts = annotate_heatmap(im, data, valfmt = valfmt, textcolors = textcolors, threshold = text_color_threshold, fontsize = fontsize_blot )
    ax.set_aspect( aspect ) # set aspect of image
    if title is not None :
        fig.suptitle( title )
    fig.tight_layout( )
    return fig, ax # return subplot axis and figure


# ### Condition Comparison Plotting

# ### Bokeh interactive plotting

# In[ ]:


basic_bokeh_label_kw = dict( x_offset = 5, y_offset = 5, text_font_size = '10pt', render_mode = 'canvas' )
def Bokeh_scatter( df, x_axis, y_axis, x_axis_type = 'linear', y_axis_type = 'linear', color = None, color_default = '#2326EE', color_log_transform = False, color_n_std_for_outliers = None, cmap = 'viridis_r', colormap_diverging = False, 
                  size = None, size_factor = 30, size_log_transform = True, size_inverse = True, size_n_std_for_outliers = 2, 
                  annotation = None, df_filtered_for_annotation = None, show_annotation_kw = basic_bokeh_label_kw, 
                  title = 'Scatter Plot', save_html = False, plot_size = 700, URL = 'https://www.google.com/search?q=@', alpha_fill_quadrants = 0.2, l_column_for_hover_tool = None, graph_folder = None ) :
    '''  General Function to draw interactive Bokeh plot  '''
    if color is not None and 'p_value' in color : # optimize setting if 'p_value' was given as one of the columns
        color_log_transform = True
    if size is not None and 'p_value' in size :
        size_log_transform, size_inverse, size_n_std_for_outliers = True, True, None
    df = deepcopy( df ).reset_index( ) # reset index so that index can be also used to plot graphs
    dict_result = df.to_dict( orient = 'list' ) # build bokeh source data dictionary from a dataframe 
    dict_result[ 'size' ] = list( np.full( len( df ), size_factor / 3 ) ) if size is None else BOKEH_transform_arr_values_into_a_range_from_0_to_1( dict_result[ size ], n_std_for_outliers = size_n_std_for_outliers, log_transform = size_log_transform, inverse = size_inverse ) * size_factor # set data for sizes of markers
    dict_result[ 'color' ] = list( np.full( len( df ), color_default ) ) if color is None else COLORMAP_GET_list_of_colors_from_list_of_values( dict_result[ color ], n_std_for_outliers = color_n_std_for_outliers, cmap = cmap, return_RGB_hex_str = True, log_transform = color_log_transform, diverging_colormap = colormap_diverging ) # set data for colors of markers
    # create source data
    source = ColumnDataSource( data = dict_result )
    # set output type (save figure as a html file)
    if save_html :
        reset_output( ) # reset output before plotting
        output_file( graph_folder + To_window_path_compatible_str( title ) + '.html' )
    else :
        output_notebook( )
    # set hovering tool according to the type of graph
    TOOLTIPS = [ ( '{x}'.format( x = x_axis ), '@{x}'.format( x = x_axis ) ), ( '{y}'.format( y = y_axis ), '@{y}'.format( y = y_axis ) ) ]
    if size is not None : TOOLTIPS.append( ( 'Size ({size})'.format( size = size ), '@{size}'.format( size = size ) ) )
    if color is not None : TOOLTIPS.append( ( 'Color ({color})'.format( color = color ), '@{color}'.format( color = color ) ) )
    if annotation is not None : TOOLTIPS.append( ( 'Annotation ({annotation})'.format( annotation = annotation ), '@{annotation}'.format( annotation = annotation ) ) )
    if l_column_for_hover_tool is not None : # if 'l_column_for_hover_tool' is given, add give list of columns to TOOLTIPs for hover tool visualization
        for col in l_column_for_hover_tool : TOOLTIPS.append( ( col, '@' + col ) )
#     title += '   color ({}) and size ({})'.format( color, size ) # add color and size annotation to the title
    p1 = figure(  tools = "box_select,lasso_select,tap,pan,box_zoom,wheel_zoom,reset,save,hover", tooltips = TOOLTIPS, title = title, x_axis_label = x_axis, y_axis_label = y_axis, x_axis_type = x_axis_type, y_axis_type = y_axis_type, plot_width = plot_size, plot_height = plot_size )
    p1.add_layout( BoxAnnotation( right = 0, bottom = 0, fill_alpha = alpha_fill_quadrants, fill_color = 'black' ) ) # add box annotations to mark boundaries
    p1.add_layout( BoxAnnotation( left = 0, top = 0, fill_alpha = alpha_fill_quadrants, fill_color = 'black' ) )
    p1.circle( x_axis, y_axis, source = source, line_color = 'black', fill_color = 'color', size = 'size' ) # plot two groups with colors
    if annotation is not None :
        taptool = p1.select( type = TapTool ) # open NCBI page with Entrez Gene ID when clicked
        taptool.callback = OpenURL( url = URL.replace( '@', '@' + annotation ) ) # call back using column given as 'annotation' column
        if df_filtered_for_annotation is not None : # if both annotation and df_filtered_for_annotation were given, show annotations on the scatter plot
            source_filtered = ColumnDataSource( data = df_filtered_for_annotation.to_dict( orient = 'list' ) ) # build bokeh source data dictionary of a dataframe after filtering  
            labels = LabelSet( x = x_axis, y = y_axis, text = annotation, level = 'glyph', source = source_filtered, ** show_annotation_kw )
            p1.add_layout( labels )
    show( p1, notebook_url = remote_jupyter_proxy_url  )


# ##### Advanced Bokeh Plotting (with Entrez Gene ID)

# In[ ]:


basic_filter = dict( thres_abs_log2fc = 0.5, thres_p_value = 0.05 )
basic_bokeh_label_kw = dict( x_offset = 5, y_offset = 5, text_font_size = '10pt', render_mode = 'canvas' )
def BOKEH_PLOT_interactive_scatter_plot_log2fc( df, x_axis = None, y_axis = None, color = 'Log2_Fold_Change', color_log_transform = False, color_n_std_for_outliers = None, size = 'p_value', size_log_transform = True, size_inverse = True, size_n_std_for_outliers = 2, x_axis_type = 'linear', y_axis_type = 'linear', show_symbol = True, show_symbol_filter_kw = basic_filter, show_symbol_kw = basic_bokeh_label_kw, title = 'Scatter Plot', save_html = False, plot_size = 700, tap_callback = 'google', size_factor = 30, color_map = 'viridis_r', colormap_diverging = False, graph_folder = None ) :
    ''' Draw interactibe Bokeh volcano plot of a given result dataframe from 'Calculate_Log2FC_p_value__A_vs_B' that calculate Log2FC
    The dataframe should have Gene_ID as its index, and all Gene_IDs has to be valid (that is, it should exist in dict_ID_Symbol)
    color of circle will be set by data of data_label given by 'color', and size of circle will be set by data of data_label given by 'size'
    'show_symbol_filter_kw' keyworded arguments for filtering entries for displaying Approved_Symboles using 'RESULT_filter'. 
    column labels of x_axis, y_axis are given thorugh arguments.    axis type, 'log' or 'linear', is also given thorugh arguments.    tap_callback : 'google' search  '''
    if x_axis is None : # if x_axis and y_axis values were not given, automatically set x and y axis values
        x_axis = df.columns.values[ 5 ]
    if y_axis is None :
        y_axis = df.columns.values[ 6 ]
    if 'p_value' in color :
        color_log_transform = True
    if 'p_value' in size :
        size_log_transform, size_inverse, size_n_std_for_outliers = True, True, None
    df = deepcopy( df )
    df = PD_Add_gene_annotation( df )
    df[ 'Gene_ID' ] = list( df.index.values )
    dict_result = df.to_dict( orient = 'list' ) # build bokeh source data dictionary from a dataframe 
    # dict_result[ 'Gene_Set_Name' ] = list( df.index.values ) # retrive Gene_Set_names
    dict_result[ 'size' ] = BOKEH_transform_arr_values_into_a_range_from_0_to_1( dict_result[ size ], n_std_for_outliers = size_n_std_for_outliers, log_transform = size_log_transform, inverse = size_inverse ) * size_factor # set data for sizes of markers
    dict_result[ 'color' ] = COLORMAP_GET_list_of_colors_from_list_of_values( dict_result[ color ], n_std_for_outliers = color_n_std_for_outliers, cmap = color_map, return_RGB_hex_str = True, log_transform = color_log_transform, diverging_colormap = colormap_diverging ) # set data for colors of markers
    # create source data
    source = ColumnDataSource( data = dict_result )
    df_filtered = RESULT_filter( df, **show_symbol_filter_kw ) 
    source_filtered = ColumnDataSource( data = df_filtered.to_dict( orient = 'list' ) ) # build bokeh source data dictionary of a dataframe after filtering  
    # set output type (save figure as a html file)
    if save_html :
        reset_output( ) # reset output before plotting
        output_file( graph_folder + To_window_path_compatible_str( title ) + '.html' )
    else :
        output_notebook( )
    # set hovering tool according to the type of graph
    TOOLTIPS = [ ( 'Gene Name (Symbol)', '@Approved_Name' ), ( '{x}'.format( x = x_axis ), '@{x}'.format( x = x_axis ) ), ( '{y}'.format( y = y_axis ), '@{y}'.format( y = y_axis ) ), ( '{size}'.format( size = size ), '@{size}'.format( size = size ) ), ( '{color}'.format( color = color ), '@{color}'.format( color = color ) ) ]
    title += '   color ({}) and size ({})'.format( color, size ) # add color and size annotation to the title
    p1 = figure(  tools = "box_select,lasso_select,tap,pan,box_zoom,wheel_zoom,reset,save,hover", tooltips = TOOLTIPS, title = title, x_axis_label = x_axis, y_axis_label = y_axis, x_axis_type = x_axis_type, y_axis_type = y_axis_type, plot_width = plot_size, plot_height = plot_size )
    fill_alpha = 0.3
    p1.add_layout( BoxAnnotation( right = 0, bottom = 0, fill_alpha = fill_alpha, fill_color = 'black' ) ) # add box annotations to mark boundaries
    p1.add_layout( BoxAnnotation( right = 0, top = 0, fill_alpha = fill_alpha, fill_color = 'black' ) )
    p1.add_layout( BoxAnnotation( left = 0, top = 0, fill_alpha = fill_alpha, fill_color = 'black' ) )
    # plot two groups with colors
    p1.circle( x_axis, y_axis, source = source, line_color = 'black', fill_color = 'color', size = 'size' )
    # open NCBI page with Entrez Gene ID when clicked
    taptool = p1.select( type = TapTool )
    if tap_callback.lower( ) == 'google' : # set URL that will be directed if marker is clicked
        URL = 'https://www.google.com/search?q=@Approved_Symbol gene wiki'
    else : # default: ncbi search
        URL = 'https://www.ncbi.nlm.nih.gov/gene/?term=@Approved_Symbol'
    taptool.callback = OpenURL( url = URL )
    # add labels showing indices on the graph of show_symbol = True
    if show_symbol :
        labels = LabelSet( x = x_axis, y = y_axis, text = 'Approved_Symbol', level = 'glyph', source = source_filtered, **show_symbol_kw )
        p1.add_layout( labels )
    show( p1 )


# In[ ]:


def BOKEH_PLOT_interactive_volcano_plot( df, x_axis = 'Log2_Fold_Change', y_axis = 'adjusted_p_value', x_axis_type = 'linear', y_axis_type = 'log', show_symbol = True, show_symbol_filter_kw = basic_filter, show_symbol_gene_ids = None, show_symbol_kw = basic_bokeh_label_kw, title = None, save_html = False, plot_size = 700, tap_callback = 'google', graph_folder = None ) :
    ''' Draw interactibe Bokeh volcano plot of a given result dataframe from 'Calculate_Log2FC_p_value__A_vs_B' that calculate Log2FC
    The dataframe should have Gene_ID as its index, and all Gene_IDs has to be valid (that is, it should exist in dict_ID_Symbol)
    column labels of x_axis, y_axis are given thorugh arguments.    axis type, log or linear, is also given thorugh arguments
    tap_callback : 'google' or 'ncbi' search of gene if clicked (default ncbi). 
    Also, annotate a given list of genes in 'show_symbol_gene_ids' if given (dafault is all genes) according to a filter given by 'show_symbol_filter_kw' '''
    dict_shorthand = dict( Log2_Fold_Change = 'Log2FC', adjusted_p_value = 'adj p', p_value = 'p val', average_all = 'avg_all' ) # define dictionary of shorthand annotation of possible x and y axes 
    df = deepcopy( df )
    df = PD_Add_gene_annotation( df )
    df[ 'Gene_ID' ] = list( df.index.values )
    dict_result = df.to_dict( orient = 'list' ) # build bokeh source data dictionary from a dataframe 
    source = ColumnDataSource( data = dict_result ) # create source data
    df_filtered = RESULT_filter( df, **show_symbol_filter_kw ) if show_symbol_gene_ids is None else RESULT_filter( PD_Subset( df, show_symbol_gene_ids ), **show_symbol_filter_kw )
    source_filtered = ColumnDataSource( data = df_filtered.to_dict( orient = 'list' ) ) # build bokeh source data dictionary of a dataframe after filtering  
    col_a, col_b = list( col for col in df.columns if 'average' in col and col != 'average_all' ) # retrive column names for condition A and B
    title = "({}) VolcanoPlot {} vs {}".format( dataset_name, col_a.replace( '_average', '' ), col_b.replace( '_average', '' ) ) if title is None else title # set title based on the names of columns
    if save_html : # set output type (save figure as a html file)
        reset_output( ) # reset output before plotting
        output_file( graph_folder + To_window_path_compatible_str( title ) + '.html' )
    else :
        output_notebook( )
    # set hovering tool according to the type of graph
    TOOLTIPS = [( '{x_axis_short}, {y_axis_short}, p val'.format( x_axis_short = dict_shorthand[ x_axis ], y_axis_short = dict_shorthand[ y_axis ] ), '@{x_axis}, @{y_axis} @p_value'.format( x_axis = x_axis, y_axis = y_axis ) ), ( 'Gene Name (Symbol)', '@Approved_Name' ), ( 'Gene_ID', '@Gene_ID' ), ( 'Cond. A vs B avg.', '@{} vs @{}'.format( col_a, col_b ) ) ]
    p1 = figure(  tools = "box_select,lasso_select,tap,pan,box_zoom,wheel_zoom,reset,save,hover", tooltips = TOOLTIPS, title = title, x_axis_label = x_axis, y_axis_label = y_axis, x_axis_type = x_axis_type, y_axis_type = y_axis_type, plot_width = plot_size, plot_height = plot_size )
    p1.circle( x_axis, y_axis, size = 8, source = source, color = 'navy', alpha = 0.3 ) # plot two groups with colors
    p1.circle( x_axis, y_axis, size = 10, source = source_filtered, color = 'navy', alpha = 1 )
    decreased_box = BoxAnnotation( right = 0, fill_alpha = 0.1, fill_color = 'green') # add box annotations # annotate upregulated genes on green background
    increased_box = BoxAnnotation( left = 0, fill_alpha = 0.1, fill_color = 'red' ) # annotate upregulated genes on red background
    p1.add_layout( decreased_box )
    p1.add_layout( increased_box )
    taptool = p1.select( type = TapTool ) # open NCBI page or Google search result with Entrez Gene ID when clicked
    if tap_callback.lower( ) == 'google' : # set URL that will be directed if marker is clicked
        URL = 'https://www.google.com/search?q=@Approved_Symbol Gene'
    elif tap_callback.lower( ) == 'gepia' :
        URL = 'http://gepia.cancer-pku.cn/detail.php?gene=@Approved_Symbol'
    else : # catch-all default: ncbi search
        URL = 'https://www.ncbi.nlm.nih.gov/gene/?term=@Gene_ID'
    taptool.callback = OpenURL( url = URL )
    if show_symbol : # add labels showing indices on the graph of show_symbol = True
        labels = LabelSet( x = x_axis, y = y_axis, text = 'Approved_Symbol', level = 'glyph', source = source_filtered, **show_symbol_kw )
        p1.add_layout( labels )
    show( p1 )


# In[ ]:


basic_filter_geneset = dict( thres_abs_log2fc = 0.5, thres_p_value = 0.05, label_p_value = 'adjusted_p_value' )
def BOKEH_PLOT_interactive_volcano_plot__Gene_Sets( df, x_axis = 'Log2_Fold_Change', y_axis = 'adjusted_p_value', x_axis_type = 'linear', y_axis_type = 'log', color = 'Log2_Fold_Change_std', color_n_std_for_outliers = 2, color_log_transform = False, colormap_diverging = False, size = 'number_valid_genes', size_log_transform = False, size_inverse = False, size_n_std_for_outliers = 2, show_symbol = True, show_symbol_filter_kw = basic_filter, title = 'Volcano Plot', save_html = False, plot_size = 700, tap_callback = 'google', size_factor = 5, color_map = 'viridis_r', graph_folder = None ) :
    '''
    Draw interactibe Bokeh volcano plot of a given result dataframe from 'Calculate_Log2FC_p_value__A_vs_B' that calculate Log2FC
    The dataframe should have Gene_ID as its index, and all Gene_IDs has to be valid (that is, it should exist in dict_ID_Symbol)
    column labels of x_axis, y_axis are given thorugh arguments.    axis type, 'log' or 'linear', is also given thorugh arguments.    tap_callback : 'google' search
    size of a circle is proportional to 'number_valid_genes' (divided by 'size_factor'), and color is related with 'coondition_A_B_std_avg', average of std of condition A samples and std of condition B samples, using 'color_map' with n_std_outliers = 2. 
    '''
    df = df.dropna( )
    dict_shorthand = dict( Log2_Fold_Change = 'Log2FC', adjusted_p_value = 'adj p', p_value = 'p val' ) # define dictionary of shorthand annotation of possible x and y axes 
    dict_result = df.to_dict( orient = 'list' ) # build bokeh source data dictionary from a dataframe 
    dict_result[ 'Gene_Set_Name' ] = list( df.index.values ) # retrive Gene_Set_names
    dict_result[ 'size' ] = BOKEH_transform_arr_values_into_a_range_from_0_to_1( dict_result[ size ], n_std_for_outliers = size_n_std_for_outliers, log_transform = size_log_transform, inverse = size_inverse ) * size_factor # set data for sizes of markers
    dict_result[ 'color' ] = COLORMAP_GET_list_of_colors_from_list_of_values( dict_result[ color ], n_std_for_outliers = color_n_std_for_outliers, color_map = color_map, return_RGB_hex_str = True, log_transform = color_log_transform, diverging_colormap = colormap_diverging ) # set data for colors of markers
    # create source data
    source = ColumnDataSource( data = dict_result )
    df_filtered = RESULT_filter( df, **show_symbol_filter_kw ) # build bokeh source data dictionary of a dataframe after filtering  
    dict_result_filtered = df_filtered.to_dict( orient = 'list' )
    dict_result_filtered[ 'Gene_Set_Name' ] = list( df_filtered.index.values ) # retrive Gene_Set_names
    source_filtered = ColumnDataSource( data = dict_result_filtered ) 
    # set output type (save figure as a html file)
    if save_html :
        reset_output( ) # reset output before plotting
        output_file( graph_folder + To_window_path_compatible_str( title ) + '.html' )
    else :
        output_notebook( )
    # set hovering tool according to the type of graph
    TOOLTIPS = [( '{x_axis_short}, {y_axis_short}'.format( x_axis_short = dict_shorthand[ x_axis ], y_axis_short = dict_shorthand[ y_axis ] ), '@{x_axis}, @{y_axis}'.format( x_axis = x_axis, y_axis = y_axis ) ), ( 'Cond. A vs B avg.', '@condition_A_average vs @condition_B_average' ), ( 'Name', '@Gene_Set_Name' ), ( 'n_genes', '@number_valid_genes' ) ]
    p1 = figure(  tools = "box_select,lasso_select,tap,pan,box_zoom,wheel_zoom,reset,save,hover", tooltips = TOOLTIPS, title = title, x_axis_label = x_axis, y_axis_label = y_axis, x_axis_type = x_axis_type, y_axis_type = y_axis_type, plot_width = plot_size, plot_height = plot_size )
    # plot two groups with colors
    p1.circle( x_axis, y_axis, source = source, line_color = 'black', fill_color = 'color', size = 'size' )
    # add box annotations
    decreased_box = BoxAnnotation( right = 0, fill_alpha = 0.1, fill_color = 'green') # annotate upregulated genes on green background
    increased_box = BoxAnnotation( left = 0, fill_alpha = 0.1, fill_color = 'red' ) # annotate upregulated genes on red background
    p1.add_layout( decreased_box )
    p1.add_layout( increased_box )
    # open NCBI page with Entrez Gene ID when clicked
    taptool = p1.select( type = TapTool )
    if tap_callback.lower( ) == 'google' : # set URL that will be directed if marker is clicked
        URL = 'https://www.google.com/search?q=@Gene_Set_Name pathway'
    else : # default: ncbi search
        URL = 'https://www.ncbi.nlm.nih.gov/gene/?term=@Gene_Set_Name'
    taptool.callback = OpenURL( url = URL )
    # add labels showing indices on the graph of show_symbol = True
    if show_symbol :
        labels = LabelSet( x = x_axis, y = y_axis, text = 'Gene_Set_Name', level = 'glyph', x_offset = 5, y_offset = 5, source = source_filtered, render_mode = 'canvas' )
        p1.add_layout( labels )
    show( p1 )


# In[ ]:


def BOKEH_PLOT_Scatter_Annotation( arr_coordinates, arr_labels, color_map = 'gist_rainbow', fill_alpha = 0.3, line_alpha = 0.3, size = 10, figsize = ( 1000, 700 ), title = 'plot', x_label = 'x_axis', y_label = 'y_axis', dict_scatter = dict( ), graph_folder = None, ** dict_attributes ) :
    '''  'arr_coordinates' should be numpy array with 2 columns, x and y coordinates, and 'arr_labels' should be a list-like object containing labels. Additional information can be given
    through keyworded arguments, and will be displayed via Bokeh's hover tool. bokeh circle function's additional keywords can be given through 'dict_scatter' '''
    arr_labels = arr_labels if isinstance( arr_labels, ( np.ndarray ) ) else np.array( arr_labels, dtype = object ) # convert arr_labels to a numpy object
    s_label = LIST_COUNT( arr_labels, duplicate_filter = None ).sort_values( ascending = False )
    save_html = False
    l_colors = [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ] if color_map is None else COLORMAP_GET_list_of_colors_from_list_of_values( np.arange( len( s_label ) ), return_RGB_hex_str = True, color_map = color_map ) # set color map for plotting 
    arr_colors = np.zeros_like( arr_labels )
    for label, color in zip( s_label.index, l_colors ) : # prepare an array of colors matching labels
        arr_colors[ arr_labels == label ] = color
    l_columns = [ 'x', 'y', 'label', 'color' ]
    arr_data = np.vstack( ( arr_coordinates.T, arr_labels, arr_colors ) ) # stack all data_values (coordinates, color, and additional attributes) into one array to make a ColumnDataSource 
    for col_label, data_values in dict_attributes.items( ) :
        data_values = np.array( data_values ) 
        arr_data = np.vstack( ( arr_data, data_values ) )
        l_columns.append( col_label )
    source = ColumnDataSource( data = dict( ( col, data ) for col, data in zip( l_columns, arr_data ) ) )
    if save_html :
        reset_output( ) # reset output before plotting
        output_file( graph_folder + To_window_path_compatible_str( title ) + '.html' )
    else :
        output_notebook( )
    TOOLTIPS = [ ( 'x', '@x' ), ( 'y', '@y' ), ( 'label', '@label' ) ] + list( ( "{}".format( col_label ), "@{}".format( col_label ) ) for col_label in l_columns[ 4 : ] )
    p = figure( tools = "box_select,lasso_select,tap,pan,box_zoom,wheel_zoom,reset,save,hover", tooltips = TOOLTIPS, title = title, x_axis_label = x_label, y_axis_label = y_label, x_axis_type = 'linear', y_axis_type = 'linear', plot_width = figsize[ 0 ], plot_height = figsize[ 1 ] )
    p.circle( 'x', 'y', source = source, line_color = 'black', fill_color = 'color', size = size, fill_alpha = fill_alpha, line_alpha = line_alpha, ** dict_scatter )
    show( p )


# In[ ]:


def subplot_two_genes( ax, gene_id_1, Gene_1_Symbol, gene_id_2, Gene_2_Symbol, df, df_name, normal_type = True, tumor_type = True, grade_color = False, alpha = 0.5, n_std_for_outliers = 3, add_regression_line_normal = True, regression_type = 'intercept zero', show_correl = 'N', external_gene_1_data = None ) :
    '''
    Receive a subplot and draw a two-gene plot as a subplot
    Tumor and normal samples are indicated by blue and red colors, respectively. 
    normal_type = True, tumor_type = True : when each type is set to True, the group appears on the plot
    default df_data is df_proteome
    set data_phospho to True if entries are indices of df_phosphoproteome
    if 'two_plots' = True and both 'df' and 'df_2' are given, draw two plots using both data as subplots to allow comparison  
    draw a regression line that go through an origin (with zero intercept) for normal samples if 'add_regression_line_normal' is True
    If 'regression_intercept_zero', draw an regression line that go through an origin
    'show_correl' : if 'N' or 'normal', show correl. for normal samples, 'T' or 'tumor' for tumor, 'NT' or 'TN' or 'A' for all samples.
                    if 'show_correl' = None, do not show correlation coefficient
    'regression_type' : set type of linear regression line. 'proj' for projection-based least square linear regression, 'th' for Theil-Sen regression, and 'ls' for normal least square linear regression, 0 for normal least linear regression passing the origin '''
    if gene_id_1 is None : # retrive data (pandas  Series) that will be used in plotting # if external data has been given, used the data
        grade_color = False # currently grade_color = True option is not available
        series_gene_1_data, series_gene_2_data = external_gene_1_data, df[ All_samples ].loc[ gene_id_2 ]
    else :
        series_gene_1_data, series_gene_2_data = df[ All_samples ].loc[ gene_id_1 ], df[ All_samples ].loc[ gene_id_2 ]
    # perform correlation analysis 
    if show_correl is not None : # set sample_set according to a given 'show_correl' argument
        if show_correl in dict_argument__samples_name :
            sample_set, name_samples = dict_type_sample[ dict_argument__samples_name[ show_correl ][ 0 ] ], dict_argument__samples_name[ show_correl ][ 1 ] # set sample_set and name_samples
        else : # by default if other than None is given, sample_set is normal and tumor samples 
            sample_set, name_samples = dict_type_sample[ 'normal_or_tumor' ], 'All Samples'
        gene_1_data ,gene_2_data = GET_non_NaN_without_outliers( series_gene_1_data[ sample_set ].values, series_gene_2_data[ sample_set ].values, n_std_for_outliers = n_std_for_outliers ) # retrive data for two genes and remove outliers and NaN values
        correl_p = stats.pearsonr( gene_1_data, gene_2_data )[ 0 ] # calculate Pearson's correlation
        correl_s = stats.spearmanr( gene_1_data, gene_2_data ).correlation # calculate Spearman's correlation
    # retrive data and plot the samples of each type
    if normal_type :
        gene_1_normal = series_gene_1_data[ N_samples ].values # retrive data for normal samples
        gene_2_normal = series_gene_2_data[ N_samples ].values
        if add_regression_line_normal : # if 'add_regression_line_normal' is set to True, draw a regression line for normal samples
            label = '(normal samples)'
            gene_1_data, gene_2_data = GET_non_NaN_without_outliers( gene_1_normal, gene_2_normal, n_std_for_outliers = n_std_for_outliers ) # remove outliers and NaN values
            x_range = np.arange( 0, gene_1_normal.max( ), .05 ) # define x_range of a regression line
            if regression_type == 0 :
                slope = np.linalg.lstsq( np.array( [ gene_1_data ] ).T, gene_2_data, rcond = None )[ 0 ][ 0 ] # solving an equation a x = b to find the best slop for the data
                line = slope * x_range # calculate a line from x_range
                label += ' Least Square Through Origin'
            elif regression_type == 'th' :
                slope, intercept, lower, higher = stats.theilslopes( gene_2_data, gene_1_data ) # calculate theilslopes and plot regression line
                line, line_lower, line_higher = ( slope * x_range ) + intercept, ( lower * x_range ) + intercept, ( higher * x_range ) + intercept # calculate lines from x_range
                label += ' Theil-Sen'
                ax.plot( x_range, line_lower, '--', color = 'brown', alpha = 0.7 ) # draw a regression of lower and upper estimation of slopes
                ax.plot( x_range, line_higher, '--', color = 'brown', alpha = 0.7 ) 
            elif regression_type == 'ls' :
                slope, intercept, r_val, p_val, stderr = stats.linregress( gene_1_data, gene_2_data ) # perform linear regression to calculate slope and intercept
                line = ( slope * x_range ) + intercept # calculate a line from x_range
                label += ' Conventional Least Square'
            else :
                slope, intercept = LINREGRESS_Projection_Based_Least_Square_Linear_Regression( gene_1_data, gene_2_data, np.vstack( ( gene_1_data, gene_2_data ) ) ) # perform projection-based linear regression to calculate slope and intercept
                line = ( slope * x_range ) + intercept # calculate a line from x_range
                label += '\n Least Square' # Proj-based
            ax.plot( x_range, line, 'k', label = label, alpha = 0.7 ) # draw a regression
        ax.plot( gene_1_normal, gene_2_normal, 'o' , color = 'b', linewidth = 10, label = 'normal', alpha = alpha ) # plot normal samples after regression line (if it has been drawn)
    if tumor_type : 
        if grade_color : 
            categories = [ 'G1', 'G2', 'G3', 'G4' ] # set categories
            colors = [ 'burlywood', 'Coral', 'red', 'magenta' ]
            data_gene_1 = get_data__tumor_grades( gene_id_1, df ) # retrive tumor_grade-wise data
            data_gene_2 = get_data__tumor_grades( gene_id_2, df )
            for category, color in zip( categories, colors ) :
                ax.plot( data_gene_1[ category ], data_gene_2[ category ], 'o' , color = color, linewidth = 10, label = category, alpha = alpha )
        else :
            gene_1_tumor = series_gene_1_data[ T_samples ].values # plot tumor data
            gene_2_tumor = series_gene_2_data[ T_samples ].values
            ax.plot( gene_1_tumor, gene_2_tumor, 'o' , color = 'r', linewidth = 10, label = 'tumor', alpha = alpha )
    ax.set_xlabel( Gene_1_Symbol + ' protein amount' )
    ax.set_ylabel( Gene_2_Symbol + ' protein amount' )
    title = df_name + ' Data'
    if show_correl is not None : # if 'show_correl' is not None, show correlation coefficient data
        title += '\n({name}, Correl. Pearson = {c_p}, Spearman = {c_s})'.format( name = name_samples, c_p = round( correl_p, 3 ), c_s = round( correl_s, 3 ) )
    ax.set_title( title )
    ax.legend( ) # add legends to the plots
    sample_set = tumor_type * T_samples + normal_type * N_samples # define sample_set that is currently shown
    gene_1_data, gene_2_data = GET_non_NaN_without_outliers( series_gene_1_data[ sample_set ].values, series_gene_2_data[ sample_set ].values, n_std_for_outliers = n_std_for_outliers ) # retrive data that is shown and remove outliers and NaN values
    x_lim = gene_1_data.mean( ) + n_std_for_outliers * gene_1_data.std( ) # set x_lim and y_lim by using mean and std of data
    y_lim = gene_2_data.mean( ) + n_std_for_outliers * gene_2_data.std( )
    return ax, x_lim, y_lim # return a subplot and calculated x_lim and y_lim


# In[ ]:


def Plot_two_genes__normal_vs_tumor( gene_1 = 'CRBN', gene_2 = 'ACAT1', data_phospho = False, switch_genes = False, df = None, df_2 = None, df_name = 'df_1', df_2_name = 'df_2', normal_type = True, tumor_type = True, grade_color = False, alpha = 0.5, save_fig = False, n_std_for_outliers = 3, add_regression_line_normal_1 = True, add_regression_line_normal_2 = False, regression_type = 'proj', show_correl_1 = 'N', show_correl_2 = None, graph_folder = None ) :
    '''
    Create scatter plot of two genes. 
    Tumor and normal samples are indicated by blue and red colors, respectively. 
    normal_type = True, tumor_type = True : when each type is set to True, the group appears on the plot
    default df_data is df_proteome
    set data_phospho to True if entries are indices of df_phosphoproteome
    if df_2' is given, draw two plots using both 'df' and 'df_2' as subplots to allow comparison between two data
    'show_correl' : if 'N' or 'normal', show correl. for normal samples, 'T' or 'tumor' for tumor, 'NT' or 'TN' or 'A' for all samples. if 'show_correl' = None, do not show correlation
    if 'switch_genes' is True, another subplots with swithced genes will be drawn
    'regression_type' : set type of linear regression line. 'proj' for projection-based least square linear regression, 'th' for Theil-Sen regression, and 'ls' for normal least square linear regression '''
    if df is None :
        df = df_proteome
        df_name = 'Proteome'
    if type( gene_1 ) is pd.Series :
        gene_id_1, gene_id_2 = None, Gene_2_Gene_ID( gene_2 )
        Gene_1_Symbol, Gene_2_Symbol = gene_1.name, dict_ID_Symbol_simple[ gene_id_2 ]
    else :
        if data_phospho :
            df = df_phosphoproteome
            df_name = 'Phosphoproteome'
            gene_id_1, gene_id_2 = gene_1, gene_2
            Gene_1_Symbol, Gene_2_Symbol = df_phosphoproteome.loc[ gene_1 ].Gene_Symbol, df_phosphoproteome.loc[ gene_2 ].Gene_Symbol
        else :
            gene_id_1, gene_id_2 = Gene_2_Gene_ID( gene_1 ), Gene_2_Gene_ID( gene_2 )
            if gene_id_1 not in df.index.values or gene_id_2 not in df.index.values : # if gene_id_1 or gene_id_2 do not exist in df, return an error value
                return -1 
            Gene_1_Symbol, Gene_2_Symbol = dict_ID_Symbol_simple[ gene_id_1 ], dict_ID_Symbol_simple[ gene_id_2 ]
    df_name = '{dataset_name} {data_name}'.format( dataset_name = dataset_name, data_name = df_name )
    if df_2 is not None or switch_genes : # if second dataFrame has been given or genes are switched, draw two subplots
        fig, ax = plt.subplots( 1, 2, figsize = ( 15, 6 ), sharey = True, sharex = True )
        _, x_lim_1, y_lim_1 = subplot_two_genes( ax[ 0 ], gene_id_1 = gene_id_1, Gene_1_Symbol = Gene_1_Symbol, gene_id_2 = gene_id_2, Gene_2_Symbol = Gene_2_Symbol, df = df, df_name = df_name, normal_type = normal_type, tumor_type = tumor_type, grade_color = grade_color, alpha = alpha, n_std_for_outliers = n_std_for_outliers, add_regression_line_normal = add_regression_line_normal_1, regression_type = regression_type, show_correl = show_correl_1, external_gene_1_data = gene_1 )
        if df_2 is not None : 
            _, x_lim_2, y_lim_2 = subplot_two_genes( ax[ 1 ], gene_id_1 = gene_id_1, Gene_1_Symbol = Gene_1_Symbol, gene_id_2 = gene_id_2, Gene_2_Symbol = Gene_2_Symbol, df = df_2, df_name = df_2_name, normal_type = normal_type, tumor_type = tumor_type, grade_color = grade_color, alpha = alpha, n_std_for_outliers = n_std_for_outliers, add_regression_line_normal = add_regression_line_normal_2, regression_type = regression_type, show_correl = show_correl_2, external_gene_1_data = gene_1 )
        else : # if switch genes are set true, plot an another subplot with the same setting with subplot_1
            _, x_lim_2, y_lim_2 = subplot_two_genes( ax[ 1 ], gene_id_1 = gene_id_2, Gene_1_Symbol = Gene_2_Symbol, gene_id_2 = gene_id_1, Gene_2_Symbol = Gene_1_Symbol, df = df, df_name = df_name, normal_type = normal_type, tumor_type = tumor_type, grade_color = grade_color, alpha = alpha, n_std_for_outliers = n_std_for_outliers, add_regression_line_normal = add_regression_line_normal_1, regression_type = regression_type, show_correl = show_correl_1, external_gene_1_data = gene_1 )
        x_lim, y_lim = max( x_lim_1, x_lim_2 ), max( y_lim_1, y_lim_2 ) # set an overall x_lim and y_lim by selecting larger values of x_lim and y_lim
    else :
        fig, ax = plt.subplots( 1, 1, sharey = True )
        _, x_lim, y_lim = subplot_two_genes( ax, gene_id_1 = gene_id_1, Gene_1_Symbol = Gene_1_Symbol, gene_id_2 = gene_id_2, Gene_2_Symbol = Gene_2_Symbol, df = df, df_name = df_name, normal_type = normal_type, tumor_type = tumor_type, grade_color = grade_color, alpha = alpha, n_std_for_outliers = n_std_for_outliers, add_regression_line_normal = add_regression_line_normal_1, regression_type = regression_type, show_correl = show_correl_1, external_gene_1_data = gene_1 )
    #title = df_name + ' Data'
    #plt.suptitle( title )
    plt.xlim( left = 0, right = x_lim ) # set an overall x_lim and y_lim
    plt.ylim( bottom = 0, top = y_lim )
    plt.show( )
    if save_fig :
        if data_phospho :
            plt.savefig( graph_folder + To_window_path_compatible_str( df_name + '_two_gene_' + gene_1 + '_' + gene2 + '.png' ), dpi = 200 )
        else :
            plt.savefig( graph_folder + To_window_path_compatible_str( df_name + '_two_gene_' + Gene_1_Symbol + '_' + Gene_2_Symbol + '.png' ), dpi = 200 )
        plt.close( )
    return


# In[ ]:


def Change__normal_vs_tumor__Pathway_or_Gene_Set( Pathway = None, datatype = 'std' , process_type = 'average', Gene_Set = None ) :
    '''
    Calculate of ( average / sum )  of protein ( amounts / z_scores ) of a given set of genes, 
    or a set of genes of a given KEGG pathway
    Categoties are 'normal' tissue and 'tumor' tissue
    
    Parameters
    ----------------
    datatype = 'raw' or 'std' :  'raw'  =  use raw protein ratio data to calculate pathway activity  
                                 'std'  =  use Z score to calculate pathway activity 
    process_type = 'average' or 'sum' : 'average' = use an average of values of genes to assess a pathway activity
                                        'sum' = use a sum of values of genes to assess a pathway activity
    '''
    if datatype == 'std' :
        df = df_proteome_zscore
    else :
        datatype = 'raw'
        df = df_proteome
    # process input gene set or pathway
    if Pathway is not None :
        # check pathway name is valid
        if not Pathway in dict_pathway :
            print ( 'invalid pathway name' )
            return -1
        Gene_List = dict_pathway[ Pathway ][ 'genes' ]
        num_pathway_genes = len( Gene_List ) # calculate the number of genes in a pathway
    elif Gene_Set is not None :
        Gene_List = list( Gene_Set )
        num_pathway_genes =  len( Gene_List )
    else :
        print( 'invalid inputs' )
        return -1
    # select genes existing in proteome data
    existing_gene_set = Gene_List
    # calculate number of valid genes and total number of genes of the pathway or set and print the values 
    valid_gene_set = set( df_proteome[ list( dict_type_sample[ 'normal_or_tumor' ] ) ].loc[ existing_gene_set ].dropna( ).index )
    num_valid_genes = len( valid_gene_set )
    if num_valid_genes == 0 : # if there is no valid genes, end the method
        print( 'number of valid genes = 0' )
        return -1 
    coverage = round( num_valid_genes / num_pathway_genes * 100, 1 ) # calculate pathway or gene set coverage
    # data values 
    if process_type == 'average' : # compute average or sum of values to assess pathway activity
        normal_data = np.average( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
        tumor_data = np.average( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
    else : 
        normal_data = np.sum( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
        tumor_data = np.sum( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
    # perform T test
    p_value = stats.ttest_ind( normal_data, tumor_data )[ 1 ]
    if datatype == 'std' : 
        change_in_tumor = np.average( tumor_data ) - np.average( normal_data )
        change_in_tumor__type = 'z score'
    else :
        change_in_tumor = 'percent', np.average( tumor_data ) / np.average( normal_data ) * 100 - 100
        change_in_tumor__type = 'percent'
    return dict( num_valid_genes = num_valid_genes, num_pathway_genes = num_pathway_genes, datatype = datatype, p_value = p_value, 
                 change_in_tumor = change_in_tumor, change_in_tumor__type = change_in_tumor__type, process_type = process_type, 
                 coverage = coverage )


# In[ ]:


def PLOT_box_plot__heatmap__normal_vs_tumor__Pathway_or_Gene_Set( df = None, Pathway = None, datatype = 'std' , process_type = 'average', displaying = True, add_Gene_Name = False, p_val_thres = 1e-8, Gene_Set = None, Gene_Set_name = '', name_of_df = '', graph_folder = None ) :
    '''
    Draw box plot of ( average / sum )  of protein ( amounts / z_scores ) of a given set of genes, or a set of genes 
    of a given KEGG pathway
    categoties are 'normal' and 'tumor' tissues
    draw significance annotation on the plot
    Also draw heatmap of individual genes
    
    Parameters
    ----------------
    datatype = 'raw' or 'std' :  'raw'  =  use raw protein ratio data to calculate pathway activity  
                                 'std'  =  use Z score to calculate pathway activity 
    process_type = 'average' or 'sum' : 'average' = use an average of values of genes to assess a pathway activity
                                        'sum' = use a sum of values of genes to assess a pathway activity
    displaying = True : if False, does not display any text while processing
    add_Gene_Name = False : if True, add Gene_Name next to Gene_Symbol in the heatmap
    p_val_thres = 1e-10 : threshold p_value for a significant annotation on the boxplot
    Gene_Set = None : set of gene_id for gene set analysis
    Gene_Set_name = '' : name of a given gene set that will be used for plot title
    '''
    threshold_bothends = None # set parameters for heatmap plotting
    cmap = 'viridis'
    if df is None :
        if datatype == 'std' : # set data type (protein level / z_score) to be used
            df = df_proteome_zscore
            ylabel = ' of Standard Scores'
            cmap = 'bwr' # reset parameters for heatmap plotting
            threshold_bothends = 0.5
        else :
            df = df_proteome
            ylabel = ' of Protein Amounts'
    else :
        datatype == 'raw'
        ylabel = ' of Protein Amounts'
    # process input gene set or pathway
    if Pathway is not None :
        # check pathway name is valid
        if not Pathway in dict_pathway :
            print ( 'invalid pathway name' )
            return -1
        Gene_List = dict_pathway[ Pathway ][ 'genes' ]
        num_pathway_genes = len( Gene_List ) # calculate the number of genes in a pathway
        Gene_Set_name = Pathway.split( '  ', 1 )[ 1 ] # handle name of pathway that will be used in the plots 
    elif Gene_Set is not None :
        Gene_List = list( Gene_Set )
        num_pathway_genes =  len( Gene_List )
    else :
        print( 'invalid inputs' )
        return -1
    # select genes existing in proteome data
    existing_gene_set = set( df.index.values ).intersection( Gene_List )
    # calculate number of valid genes and total number of genes of the pathway and print the values 
    valid_gene_set = set( df[ list( dict_type_sample[ 'normal_or_tumor' ] ) ].loc[ existing_gene_set ].dropna( ).index )
    num_valid_genes = len( valid_gene_set )
    if num_valid_genes == 0 : # if there is no valid genes, end the method
        print( 'number of valid genes = 0' )
        return
    coverage = round( num_valid_genes / num_pathway_genes * 100, 1 ) # calculate pathway coverage
    if displaying : # display the result only when displaying = True
        print( num_valid_genes, ' / ', num_pathway_genes, ' genes exist in proteome data' )
    # Generate title of the plot
    title =  dataset_name + ' ' + name_of_df + " " + Gene_Set_name + "'\n(" + str( coverage ) + '% covered) ' + process_type +  ' of proteins (' + datatype + ')'
    ################################################ Ploting pathway-wise boxplot ####################################
    fig1, ax1 = plt.subplots( figsize = ( 6, 5 ) )
    ax1.set_title( title )
    # data values 
    if process_type == 'average' : # compute average or sum of values to assess pathway activity
        ylabel = 'Average' + ylabel # edit y_label for the boxplot
        normal_data = np.average( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
        tumor_data = np.average( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
    else : 
        ylabel = 'Sum' + ylabel # edit y_label for the boxplot
        normal_data = np.sum( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
        tumor_data = np.sum( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ valid_gene_set ].iloc[ :, 1: ].values, axis = 0 )
    # perform T test
    p_value = stats.ttest_ind( normal_data, tumor_data )[ 1 ]
    p_value_sig = format( p_value, '.3g') # number of significant digits = 3 # string
    if displaying : # display the results only when displaying = True
        print( 'p value : ' + p_value_sig ) # print p_value
        # print changed amount of pathway activity
        if datatype == 'std' : 
            print ( ylabel, np.average( tumor_data ) - np.average( normal_data ), ' changed in tumor' )
        else :
            print ( ylabel, np.average( tumor_data ) / np.average( normal_data ) * 100 - 100, '% changed in tumor' )
    # set Y range
    y_max = max( np.max( normal_data ), np.max( tumor_data ) ) 
    if datatype == 'std' :
        y_min = min( np.min( normal_data ), np.min( tumor_data ) )
    else :
        y_min = 0
    data = [ normal_data, tumor_data ]
    labels = ['Normal', 'Tumor']
    # plot the box plot
    ax1.boxplot( data, notch = True, labels = labels, meanline = True, showmeans = True )
    ax1.set_ybound( lower = y_min * 1.1 , upper = y_max * 1.1 + 0.2 )
    ax1.set_xlabel('Sample Type')
    ax1.set_ylabel( ylabel )
    # add significance annotation
    if p_value < p_val_thres :
        x1, x2 = 1, 2  # columns 
        y, h, col = y_max * 1.05, y_max * 0.02, 'b'
        plt.plot( [ x1, x1, x2, x2 ], [ y, y + h, y + h, y ], lw=1, c = col )
        plt.text( (x1 + x2) * .5, y + h, 'p = ' + p_value_sig , ha='center', va='bottom', color=col)
    # save and close the figure
    plt.subplots_adjust( top = 0.88 )
    plt.savefig( graph_folder + title.replace('/', '_').replace( '\n', ' ' ) + '.png' , dpi = 200 ) # save figure at graph folder
    plt.close()
    ################################################ Ploting gene-wise heatmap ####################################
    # convert the list of valid Gene_ID to list of Gene_Symbols 
    list_valid_gene = sorted( list( valid_gene_set ) ) # convert gene set to an ordered data structure, a list
    column_labels = List_Gene_ID__2__List_Gene_Symbol( list_valid_gene, add_Gene_Name )
    # set figure size depending on the number of genes, set title
    fig, ax = plt.subplots( figsize=( 8 + int( len( column_labels ) / 2.5 ), 3 + 3 * int( add_Gene_Name )  ) )
    fig.suptitle( title + ' heatmap', fontsize = 14 ) # display title name as a line 
    # Turn off grid
    ax.grid( False )
    # create data
    if process_type == 'average' : # compute average or sum of values to assess pathway activity
        normal_data = np.average( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ list_valid_gene ].iloc[ :, 1: ].values, axis = 1 )
        tumor_data = np.average( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ list_valid_gene ].iloc[ :, 1: ].values, axis = 1 )
    else : 
        normal_data = np.sum( df[ list( dict_type_sample[ 'normal' ] ) ].loc[ list_valid_gene ].iloc[ :, 1: ].values, axis = 1 )
        tumor_data = np.sum( df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ list_valid_gene ].iloc[ :, 1: ].values, axis = 1 )
    # draw heatmap
    im, cbar = heatmap( np.vstack( ( normal_data, tumor_data ) ), labels, 
                       column_labels, ax = ax, cmap = cmap, cbarlabel = ylabel )
    if len( list_valid_gene ) < 15 :
        texts = annotate_heatmap(im, valfmt="{x:.2f}", threshold_bothends = threshold_bothends )
    else :
        texts = annotate_heatmap(im, valfmt="{x:.1f}", threshold_bothends = threshold_bothends )
    fig.tight_layout( )
    # save fiqure
    fig.savefig( graph_folder + title.replace('/', '_').replace('\n', ' ') + '_heatmap.png' , dpi = 150 )
    plt.close()
    return


# In[ ]:


def PLOT_box_plot__normal_vs_tumor__Gene( Gene, datatype = 'std', displaying = True, p_val_thres = 1e-8, df = None, name_of_df = 'protein level', graph_folder = None ) :
    '''
    draw whisker plot of protein level of a given gene ID or gene_symbol
    categoties are 'normal' tissue and 'tumor' tissue
    draw significance annotation on the plot
    
    Parameters
    ----------------
    datatype = 'raw' or 'std' :  'raw'  =  use raw protein ratio data to calculate pathway activity  
                                 'std'  =  use Z score to calculate pathway activity 
    displaying = True : if False, does not display any text while processing
    p_val_thres = 1e-10 : threshold p_value for significant annotation on the boxplot
    '''
    if df is None :
        df = df_proteome
    # retrive Gene_ID and check Gene is valid
    Gene_ID = Gene_2_Gene_ID( Gene )
    if Gene_ID == -1 :
        return -1
    # set data type (protein level / z_score) to be used
    if datatype == 'std' : 
        df = df_proteome_zscore
        ylabel = 'Standard Score'
    else :
        # df = df_proteome
        ylabel = 'Protein Amount'
    # Generate title of the plot
    title =  dataset_name + "\n" + dict_ID_Symbol[ Gene_ID ][ 1 ] + ' ( ' + dict_ID_Symbol_simple[ Gene_ID ] + ' )\n ' + name_of_df + ' (' + datatype + ')'
    ################################################ Ploting pathway-wise boxplot ####################################
    fig1, ax1 = plt.subplots( figsize = ( 6, 5) )
    ax1.set_title( title )
    # data values 
    normal_data = df[ list( dict_type_sample[ 'normal' ] ) ].loc[ Gene_ID ].dropna( ).values
    tumor_data = df[ list( dict_type_sample[ 'tumor' ] ) ].loc[ Gene_ID ].dropna( ).values
    # perform T test
    p_value = stats.ttest_ind( normal_data, tumor_data )[ 1 ]
    p_value_sig = format( p_value, '.3g') # number of significant digits = 3 # string
    if displaying : # display the results only when displaying = True
        print( 'p value : ' + p_value_sig ) # print p_value
        # print changed amount of pathway activity
        if datatype == 'std' : 
            print ( ylabel, np.average( tumor_data ) - np.average( normal_data ), ' changed in tumor' )
        else :
            print ( ylabel, np.average( tumor_data ) / np.average( normal_data ) * 100 - 100, '% changed in tumor' )
    # set Y range
    y_max = max( np.max( normal_data ), np.max( tumor_data ) ) 
    if datatype == 'std' :
        y_min = min( np.min( normal_data ), np.min( tumor_data ) )
    else :
        y_min = 0
    data = [ normal_data, tumor_data ]
    labels = ['Normal', 'Tumor']
    # plot the box plot
    ax1.boxplot( data, notch = True, labels = labels, meanline = True, showmeans = True )
    ax1.set_ybound( lower = y_min * 1.1 , upper = y_max * 1.1 + 0.2 )
    ax1.set_xlabel('Sample Type')
    ax1.set_ylabel( ylabel )
    # add significance annotation
    if p_value < p_val_thres :
        x1, x2 = 1, 2  # columns 
        y, h, col = y_max * 1.05, y_max * 0.02, 'b'
        plt.plot( [ x1, x1, x2, x2 ], [ y, y + h, y + h, y ], lw=1, c = col )
        plt.text( (x1 + x2) * .5, y + h, 'p = ' + p_value_sig , ha='center', va='bottom', color=col)
    # save and close the figure
    plt.subplots_adjust( top = 0.88 )
    plt.savefig( graph_folder + title.replace( '\n', ' ' ) + '.png' , dpi = 200 ) # save figure at graph folder
    plt.close()


# In[ ]:


def PLOT_whisker_plot__tumor_grades__Gene( Gene, df = None ) :
    '''
    draw whisker plot of protein level of a given gene ID or gene_symbol
    categoties are based on normal, G1, G2, G3, G4 tumor grades
    '''
    # retrive Gene_ID and check Gene is valid
    Gene_ID = Gene_2_Gene_ID( Gene )
    if Gene_ID == -1 :
        return -1
    # Generate title of the plot
    title = dataset_name + ' ' + dict_ID_Symbol_simple[ Gene_ID ] + ' protein level with tumor progression'
    # define categories
    categories = ['normal', 'G1', 'G2', 'G3', 'G4']
    # data values 
    normal_data = df[ list( dict_tumor_grade__sample[ categories[ 0 ] ] ) ].loc[ Gene_ID ].values
    G1_data = df[ list( dict_tumor_grade__sample[ categories[ 1 ] ] ) ].loc[ Gene_ID ].values
    G2_data = df[ list( dict_tumor_grade__sample[ categories[ 2 ] ] ) ].loc[ Gene_ID ].values
    G3_data = df[ list( dict_tumor_grade__sample[ categories[ 3 ] ] ) ].loc[ Gene_ID ].values
    G4_data = df[ list( dict_tumor_grade__sample[ categories[ 4 ] ] ) ].loc[ Gene_ID ].values
    # Create dataSources
    source_normal = ColumnDataSource( data = dict( data = normal_data, label = [ 'normal' for i in range( len( normal_data ) ) ] ) )
    source_G1 = ColumnDataSource( data = dict( data = G1_data, label = [ 'G1' for i in range( len( G1_data ) ) ] ) )
    source_G2 = ColumnDataSource( data = dict( data = G2_data, label = [ 'G2' for i in range( len( G2_data ) ) ] ) )
    source_G3 = ColumnDataSource( data = dict( data = G3_data, label = [ 'G3' for i in range( len( G3_data ) ) ] ) )
    source_G4 = ColumnDataSource( data = dict( data = G4_data, label = [ 'G4' for i in range( len( G4_data ) ) ] ) )
    # Create a Bokeh plot
    output_notebook()
    p = figure(plot_width=650, x_range = categories, title = title)
    p.circle( y ='data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_normal, color = 'blue', alpha = 0.5 )
    p.circle( y ='data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_G1, color = 'red', alpha = 0.5 )
    p.circle( y ='data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_G2, color = 'red', alpha = 0.5 )
    p.circle( y ='data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_G3, color = 'red', alpha = 0.5 )
    p.circle( y ='data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_G4, color = 'red', alpha = 0.5 )
    # calculate standard deviation data
    base, lower, upper = [], [], []
    # define two sets of data for iteration
    sets_data = [ normal_data, G1_data, G2_data, G3_data, G4_data ]
    for data in sets_data:
        mean = np.average( data)
        std = np.std( data)
        lower.append( mean - std)
        upper.append( mean + std)
    # convert standard deviation data into dataSource
    source_error = ColumnDataSource( data = dict( base = categories, lower = lower, upper = upper))
    # add standard deviation data
    p.add_layout( Whisker(source = source_error, base = 'base', upper = 'upper', lower = 'lower') )
    show( p )
    return


# In[ ]:


def PLOT_whisker_plot__tumor_stages__Gene( Gene, df = None ) :
    '''
    draw whisker plot of protein level of a given gene ID or gene_symbol
    categoties are based on normal, T1a, T1b, T2a, T2b, T3, T3a, T3b, and T4  tumor grades
    '''
    # retrive Gene_ID and check Gene is valid
    Gene_ID = Gene_2_Gene_ID( Gene )
    if Gene_ID == -1 :
        return -1
    # Generate title of the plot
    title = dataset_name + ' ' + dict_ID_Symbol_simple[ Gene_ID ] + ' protein level with tumor invasions (tumor stages)'
    # define categories
    categories = [ 'normal', 'T1a', 'T1b', 'T2a', 'T2b', 'T3', 'T3a', 'T3b', 'T4'] 
    # data values 
    normal_data = df[ list( dict_tumor_stage__sample[ categories[ 0 ] ] ) ].loc[ Gene_ID ].values
    T1a_data = df[ list( dict_tumor_stage__sample[ categories[ 1 ] ] ) ].loc[ Gene_ID ].values
    T1b_data = df[ list( dict_tumor_stage__sample[ categories[ 2 ] ] ) ].loc[ Gene_ID ].values
    T2a_data = df[ list( dict_tumor_stage__sample[ categories[ 3 ] ] ) ].loc[ Gene_ID ].values
    T2b_data = df[ list( dict_tumor_stage__sample[ categories[ 4 ] ] ) ].loc[ Gene_ID ].values
    T3_data = df[ list( dict_tumor_stage__sample[ categories[ 5 ] ] ) ].loc[ Gene_ID ].values
    T3a_data = df[ list( dict_tumor_stage__sample[ categories[ 6 ] ] ) ].loc[ Gene_ID ].values
    T3b_data = df[ list( dict_tumor_stage__sample[ categories[ 7 ] ] ) ].loc[ Gene_ID ].values
    T4_data = df[ list( dict_tumor_stage__sample[ categories[ 8 ] ] ) ].loc[ Gene_ID ].values
    # Create dataSources
    source_normal = ColumnDataSource( data = dict( data = normal_data, label = [ 'normal' for i in range( len( normal_data ) ) ] ) )
    source_T1a = ColumnDataSource( data = dict( data = T1a_data, label = [ 'T1a' for i in range( len( T1a_data ) ) ] ) )
    source_T1b = ColumnDataSource( data = dict( data = T1b_data, label = [ 'T1b' for i in range( len( T1b_data ) ) ] ) )
    source_T2a = ColumnDataSource( data = dict( data = T2a_data, label = [ 'T2a' for i in range( len( T2a_data ) ) ] ) )
    source_T2b = ColumnDataSource( data = dict( data = T2b_data, label = [ 'T2b' for i in range( len( T2b_data ) ) ] ) )
    source_T3 = ColumnDataSource( data = dict( data = T3_data, label = [ 'T3' for i in range( len( T3_data ) ) ] ) )
    source_T3a = ColumnDataSource( data = dict( data = T3a_data, label = [ 'T3a' for i in range( len( T3a_data ) ) ] ) )
    source_T3b = ColumnDataSource( data = dict( data = T3b_data, label = [ 'T3b' for i in range( len( T3b_data ) ) ] ) )
    source_T4 = ColumnDataSource( data = dict( data = T4_data, label = [ 'T4' for i in range( len( T4_data ) ) ] ) )
    # Create a Bokeh plot
    output_notebook()
    p = figure(plot_width=650, x_range = categories, title = title)
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_normal, color = 'blue', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T1a, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T1b, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T2a, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T2b, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T3, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T3a, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T3b, color = 'red', alpha = 0.5 )
    p.circle( y = 'data', x = jitter( 'label', width = 0.4, range = p.x_range),  source = source_T4, color = 'red', alpha = 0.5 )
    # calculate standard deviation data
    base, lower, upper = [], [], []
    # define two sets of data for iteration
    sets_data = [ normal_data, T1a_data, T1b_data, T2a_data, T2b_data, T3_data, T3a_data, T3b_data, T4_data ]
    for data in sets_data:
        mean = np.average( data)
        std = np.std( data)
        lower.append( mean - std)
        upper.append( mean + std)
    # convert standard deviation data into dataSource
    source_error = ColumnDataSource( data = dict( base = categories, lower = lower, upper = upper))
    # add standard deviation data
    p.add_layout( Whisker(source = source_error, base = 'base', upper = 'upper', lower = 'lower') )
    show( p)
    return


# ### Functions for Matplotlib

# In[ ]:


def Matplotlib_change_box_plot_color( ax, edge_color, fill_color = None, alpha = None ):
    ''' A method that can edit color of boxplot patches in a given subplot  '''
    for element in [ 'boxes', 'whiskers', 'fliers', 'caps' ] : # 'means', 'medians',
        plt.setp( ax[ element ], color = edge_color )
    if fill_color == 'same' :
        pass
    elif fill_color is not None :
        for patch in ax[ 'boxes' ]:
            patch.set_facecolor( fill_color )  
    else :
        for patch in ax[ 'boxes' ]: # if fill_color has not been given, remove the color of box
            patch.set_fill( False ) 
    if alpha is not None :
        for element in [ 'boxes', 'whiskers', 'fliers', 'caps' ] : # 'means', 'medians',
            plt.setp( ax[ element ], alpha = alpha )


# ### Plotting functions with PCA, t-SNE, and clustering analysis 

# In[ ]:


def PLOT_PCA_plot( df, annotate_labels = None ) :
    ''' df : ( Genes = index, samples = column ). Draw 2D PCA plot of samples. Annotate samples in 'annotate_labels' if other than None value is given.  '''
    data = df.values        
    fig, ax = plt.subplots( 1, 1 )
    for_PCA_plot = PCA( n_components = 2 ) # PCA transformation to three components 
    PCA_result = for_PCA_plot.fit_transform( data.T )
    x_coord, y_coord = PCA_result.T
    ax.scatter( x_coord, y_coord )
    ax.set_xlabel( 'PCA_1' )
    ax.set_ylabel( 'PCA_2' )
    if annotate_labels is not None :
        x_step, y_step = ( x_coord.max( ) - x_coord.min( ) ) / 100, ( y_coord.max( ) - y_coord.min( ) ) / 100 # define x and y step for annotation
        mask = GET_MASK_of_intersection( df.columns.values, annotate_labels )
        for x, y, label in zip( x_coord[ mask ], y_coord[ mask ], df.columns.values[ mask ] ) :
            ax.annotate( label, ( x, y ), xytext = ( x + x_step, y + y_step ) )


# In[ ]:


def PLOT_3D_PCA_plot( data ) :
    ''' data : ( Genes = index, samples = column ) 
    Draw 3D PCA plot of samples '''
    if type( data ) is pd.DataFrame :
        data = data.values        
    fig = plt.figure( )
    ax = fig.add_subplot( 111, projection = '3d' )
    for_PCA_plot = PCA( n_components = 3 ) # PCA transformation to three components 
    PCA_result = for_PCA_plot.fit_transform( data.T )
    x_coord, y_coord, z_coord = PCA_result.T
    ax.scatter( x_coord, y_coord, z_coord )
    ax.set_xlabel('PCA_1')
    ax.set_ylabel('PCA_2')
    ax.set_zlabel('PCA_3')


# In[ ]:


def PLOT_Embedding_2D( df, arr_2D_embedding = None, n_features_pca = 100, method_embedding = 'UMAP', Kmean_clusters = None, perplexity = 50, metric = 'correlation', dbscan_eps_min_samples = ( 0.5, 7 ), plot_samples_not_clustered = True, annotate_labels = None, individual_annotation = False, figsize = ( 8, 6 ), ** dict_embedding ) :
    ''' data : ( Genes = index, samples = column ) 
    Draw 2D tSNE or UMAP transformed embedding of a given dataframe arroding to 'method_embedding'
    return 'arr_2D_embedding' and cluster_labels for further characterization '''
    data = df.values # retrive a numpy array data
    if arr_2D_embedding is None : # if previously transformed t-SNE data is not given, by using PCA and t-SNE, perform t-SNE transformation 
        data = data.T
        reducer_pca = PCA( n_components = np.min( [ n_features_pca, data.shape[ 0 ] - 1, data.shape[ 1 ] - 1 ] ) ) # n_component of PCA transformation is minimum between 'n_features_pca', number of samples, or number of features in data
        arr_pca = reducer_pca.fit_transform( data )
        if 'sne' in method_embedding.lower( ) :
            reducer = TSNE( n_components = 2, perplexity = perplexity, ** dict_embedding ) # perform t-SNE transformation with 2 components
        elif 'umap' == method_embedding.lower( ) :
            reducer = umap.UMAP( metric = metric, ** dict_embedding )
        arr_2D_embedding = reducer.fit_transform( arr_pca ) # columns = Genes, indices = Cell_IDs
    fig, ax = plt.subplots( figsize = figsize )  
    if Kmean_clusters is not None : # if parameter for Kmean cluster has been given, use Kmean to cluster TSNE result
        kmeans = KMeans( n_clusters = Kmean_clusters, random_state = 0 ).fit( arr_2D_embedding ) # Kmean cluster data using Kmean_cluster argument
        n_clusters, cluster_labels = Kmean_clusters, kmeans.labels_
    else : # by default, use DBSCAN method implemented in sklearn to estimates optimal number of clusters and cluster ths samples
        db = DBSCAN( eps = dbscan_eps_min_samples[ 0 ], min_samples = dbscan_eps_min_samples[ 1 ] ).fit( arr_2D_embedding )
        n_clusters, cluster_labels = len( set( db.labels_ ) ) - ( 1 if - 1 in db.labels_ else 0 ), db.labels_
    for index_cluster in range( -1 * plot_samples_not_clustered, n_clusters ) : # draw markers for samples in each cluster
        x_coord, y_coord = arr_2D_embedding[ cluster_labels == index_cluster ].T 
        ax.scatter( x_coord, y_coord, label = 'cluster {index_cluster}'.format( index_cluster = index_cluster ), alpha = 0.7 ) # for each cluster, plot markers with different colors
    if annotate_labels is not None : # if list of samples for annotation has been given
        samples = df.columns.values # retrive list of samples
        if individual_annotation : # if individial_annotation is set to True, annotate individual samples with unique labels 
            annotate_labels = np.array( annotate_labels )[ GET_MASK_of_intersection( annotate_labels, samples ) ] # retrive valid samples from the given list while preserving the order
            for sample in annotate_labels :
                index = np.where( samples == sample )[ 0 ][ 0 ]
                x_coord_selected, y_coord_selected = TSNE_transformed[ index ]
                ax.scatter( x_coord_selected, y_coord_selected, label = sample, linewidths = 10 ) # annotate individual samples with unique label and (possibly) color
        else : # individial_annotation is set to False, annotate all valid samples in annotate_labels as a whole
            x_coord_selected, y_coord_selected = arr_2D_embedding[ GET_MASK_of_intersection( samples, annotate_labels ) ].T
            ax.scatter( x_coord_selected, y_coord_selected, label = 'Annotated Samples', color = 'red', linewidths = 5 )
    ax.set_xlabel( '{}_1'.format( method_embedding ) )
    ax.set_ylabel( '{}_2'.format( method_embedding ) )
    plt.legend( )
    if annotate_labels is not None :
        x_step, y_step = ( x_coord.max( ) - x_coord.min( ) ) / 100, ( y_coord.max( ) - y_coord.min( ) ) / 100 # define x and y step for annotation
        mask = GET_MASK_of_intersection( df.columns.values, annotate_labels )
        for x, y, label in zip( x_coord[ mask ], y_coord[ mask ], df.columns.values[ mask ] ) :
            ax.annotate( label, ( x, y ), xytext = ( x + x_step, y + y_step ) )
    return arr_2D_embedding, cluster_labels # return TSNE-Transformed data for further optimization of cluster number
PLOT_tSNE_plot = PLOT_Embedding_2D


# In[ ]:


def PLOT_3D_tSNE_plot( df, TSNE_transformed = None, Kmean_clusters = None, dbscan_eps_min_samples = ( 0.5, 7 ), TSNE_learning_rate = 100, plot_samples_not_clustered = True, sample_list_for_annotation = None, individual_annotation = False ) :
    ''' data : ( Genes = index, samples = column ) 
    Draw 3D PCA plot of samples '''
    data = df.values # retrive a numpy array data
    if TSNE_transformed is None : # if previously transformed t-SNE data is not given, by using PCA and t-SNE, perform t-SNE transformation 
        data = data.T
        PCA_for_tSNE = PCA( n_components = np.min( [ 50, len( data ) ] ) ) # n_component of PCA transformation is minimum between 50 and number of samples in data
        PCA_transformed = PCA_for_tSNE.fit_transform( data )
        for_3D_Plot = TSNE( n_components = 3, learning_rate = TSNE_learning_rate ) # perform t-SNE transformation with 3 components
        TSNE_transformed = for_3D_Plot.fit_transform( PCA_transformed )
    fig = plt.figure( )
    ax = fig.add_subplot( 111, projection = '3d' ) 
    if Kmean_clusters is not None : # if parameter for Kmean cluster has been given, use Kmean to cluster TSNE result
        kmeans = KMeans( n_clusters = Kmean_clusters, random_state = 0 ).fit( TSNE_transformed ) # Kmean cluster data using Kmean_cluster argument
        n_clusters, cluster_labels = Kmean_clusters, kmeans.labels_
    else : # by default, use DBSCAN method implemented in sklearn to estimates optimal number of clusters and cluster ths samples
        db = DBSCAN( eps = dbscan_eps_min_samples[ 0 ], min_samples = dbscan_eps_min_samples[ 1 ] ).fit( TSNE_transformed )
        n_clusters, cluster_labels = len( set( db.labels_ ) ) - ( 1 if - 1 in db.labels_ else 0 ), db.labels_
    for index_cluster in range( -1 * plot_samples_not_clustered, n_clusters ) : # draw markers for samples in each cluster
        x_coord, y_coord, z_coord = TSNE_transformed[ cluster_labels == index_cluster ].T 
        ax.scatter( x_coord, y_coord, z_coord, label = 'cluster {index_cluster}'.format( index_cluster = index_cluster ) ) # for each cluster, plot markers with different colors
    if sample_list_for_annotation is not None : # if list of samples for annotation has been given
        samples = df.columns.values # retrive list of samples
        if individual_annotation : # if individial_annotation is set to True, annotate individual samples with unique labels 
            sample_list_for_annotation = np.array( sample_list_for_annotation )[ GET_MASK_of_intersection( sample_list_for_annotation, samples ) ] # retrive valid samples from the given list while preserving the order
            for sample in sample_list_for_annotation :
                index = np.where( samples == sample )[ 0 ][ 0 ]
                x_coord, y_coord, z_coord = TSNE_transformed[ index ]
                ax.scatter( x_coord, y_coord, z_coord, label = sample, linewidths = 10 ) # annotate individual samples with unique label and (possibly) color
        else : # individial_annotation is set to False, annotate all valid samples in sample_list_for_annotation as a whole
            x_coord, y_coord, z_coord = TSNE_transformed[ GET_MASK_of_intersection( samples, sample_list_for_annotation ) ].T
            ax.scatter( x_coord, y_coord, z_coord, label = 'Annotated Samples', color = 'k', linewidths = 10 )
    ax.set_xlabel('t-SNE_1')
    ax.set_ylabel('t-SNE_2')
    ax.set_zlabel('t-SNE_3')
    plt.legend( )
    return TSNE_transformed # return TSNE-Transformed data for further optimization of cluster number


# In[ ]:


def PLOT_DimRed_Gene( df, arr_2d, gene, linewidth = 0.7, radius = 0.5, color_map = 'viridis', figsize = ( 8, 6 ), geneset_name = 'Collection of Genes', ** dict_circle_or_plot ) :
    ''' Plot Intensity of a given Gene or the average intensity of a given list of genes for each sample at x y coordinates from a dimentionality reduction. To use plt.Circle function, put 'circle' = True
    The name of gene set can be given by the 'geneset_name' argument '''
    s_data, name_gene = ( PD_Locate_gene( df, gene ), gene ) if isinstance( gene, ( str, int, float ) ) else ( PD_Subset( df, gene ).mean( axis = 0 ), geneset_name ) # retrive a data of a gene or avarage of a given list of genes
    data = np.vstack( [ arr_2d.T.astype( object ), COLORMAP_GET_list_of_colors_from_list_of_values( s_data.values, return_RGB_hex_str = True, color_map = color_map ).astype( object ) ] )
    fig, ax = plt.subplots( 1, figsize = figsize )
    flag_plot_circle = dict_circle_or_plot.pop( 'circle', False ) # if circle = True argument was given, set 'flag_plot_circle' to True
    if flag_plot_circle :
        for dim_1, dim_2, color in data.T :
            ax.add_artist( plt.Circle( ( dim_1, dim_2 ), radius = radius, linewidth = linewidth, color = color, fill = False, ** dict_circle_or_plot ) )
    else :
        ax.scatter( data[ 0 ], data[ 1 ], linewidth = linewidth, color = data[ 2 ], ** dict_circle_or_plot )
    x_lim, y_lim = np.vstack( [ np.min( arr_2d, axis = 0 ), np.max( arr_2d, axis = 0 ) ] ).T * 1.1 # retrive x and y limits of the plot
    MATPLOTLIB_basic_configuration( x_lim = tuple( x_lim ), y_lim = tuple( y_lim ), title = name_gene )


# ### Plotting Functions for Advanced Analysis

# In[ ]:


def PLOT_NORMALIZE_boxplots_organelle_amount_along_with_Gene_data( df, Organelle_Gene_Set, Gene = None, plot_normalized = False, alpha = 0.5, box_width = 0.0035, n_skip_samples = 1, show_ref = False, Gene_Set_description = 'Mitochondrial' ) :
    ''' subset 'df' with 'Organelle_Gene_Set' to represent an organelle and RLE normalize the data
        With an estimated amount of organelle in each sample, plots boxplots of each samples, along with scatter plots using Gene data (set N_samples as pseudo_reference )
        if 'Gene' is given. If 'plot_normalized' is set to True, boxplots and the scatter plot are drawn using the normalized data 
        external data can be given through 'Gene' argument as a pandas series where index is Sample_ID '''
    df_all_input_genes = deepcopy( df )
    df_normalized, series_organelle_amount = NORMALIZE_Relative_Log_Expression( df, N_samples, Organelle_Gene_Set, return_normalized_ratios = True, scale_all_genes = True ) # retrive normalizaed data and estimated mitochondria amount from RLE
    df = PANDAS_Subset( df, Organelle_Gene_Set ).dropna( ) # set a dataframe to display
    num_Organelle_Genes = len( df ) # retrive number of 'Organelle_Gene_Set' genes that exist in df 
    selection_mask = np.array( np.zeros( len( series_organelle_amount ) ), dtype = bool )
    selection_mask[ : : n_skip_samples ] = True # mask that will select every n_skip_samples-th samples from the sorted
    list_selected_sorted_samples = series_organelle_amount.sort_values( )[ selection_mask ].index.values
    list_tumor_samples_sorted = list( sample for sample in list_selected_sorted_samples if sample in T_samples )
    list_normal_samples_sorted = list( sample for sample in list_selected_sorted_samples if sample in N_samples ) # retrive sorted list of tumor samples and normal samples
    list_ref_samples_sorted = list( sample for sample in list_selected_sorted_samples if sample in Ref_samples ) # retrive sorted list of tumor samples and normal samples
    if plot_normalized : # if 'plot_normalized' is set to True, use normalized organelle data to plot the data
        df_all_input_genes = df_normalized # set a dataframe that will be used for plotting
        df = PANDAS_Subset( df_normalized, Organelle_Gene_Set ) # subset 'Organelle_Gene_Set' fron normalized data
    tumor_data = df[ list_tumor_samples_sorted ].values # retrive data that will be used for boxplots
    normal_data = df[ list_normal_samples_sorted ].values
    ref_data = df[ list_ref_samples_sorted ].values
    if Gene is not None : # set colors for plots if Gene is given
        fill_color, ref_bp_color, ref_gene_color, N_bp_color, N_gene_color, T_bp_color, T_gene_color = None, 'chartreuse', 'g', 'cyan', 'b', 'darkorange', 'r'
    else : # if Gene is not given, set colors and boxplot_width for plots accordingly.
        box_width *= 3 
        fill_color, ref_bp_color, N_bp_color, T_bp_color = 'same', 'g', 'b', 'r'
    tumor_pos = series_organelle_amount[ list_tumor_samples_sorted ].values # retrive the positions of boxplots of tumor, normals, and ref
    normal_pos = series_organelle_amount[ list_normal_samples_sorted ].values
    ref_pos = series_organelle_amount[ list_ref_samples_sorted ].values
    fig1, ax1 = plt.subplots( )# start plotting
    if show_ref : # show boxplots of reference samples, if 'show_ref' is set to True
        bp_ref = ax1.boxplot( ref_data, notch = True, meanline = False, showmeans = False, positions = ref_pos, widths = box_width, showfliers = False, patch_artist = True, medianprops = dict( alpha = alpha, color = ref_bp_color ) )
        Matplotlib_change_box_plot_color( bp_ref, ref_bp_color, fill_color, alpha = alpha )
    bp_tumor = ax1.boxplot( tumor_data, notch = True, meanline = False, showmeans = False, positions = tumor_pos, widths = box_width, showfliers = False, patch_artist = True, medianprops = dict( alpha = alpha, color = T_bp_color ) )
    Matplotlib_change_box_plot_color( bp_tumor, T_bp_color, fill_color, alpha = alpha )
    bp_normal = ax1.boxplot( normal_data, notch = True, meanline = False, showmeans = False, positions = normal_pos, widths = box_width, showfliers = False, patch_artist = True, medianprops = dict( alpha = alpha, color = N_bp_color ) )
    Matplotlib_change_box_plot_color( bp_normal, N_bp_color, fill_color, alpha = alpha )
    is_normalized = False # set a flag for displying average-scaled gene_data 
    if Gene is not None : # if Gene is given, retrive data of gene and plot the dataalong with boxplots representing organelles
        if type( Gene ) is pd.Series : # if external data was given through 'Gene' argument, use them to plot a data 
            is_proteome_data = False # a flag that will be used for labels of plots
            data_gene = deepcopy( Gene ) # copy the given external dta
            Gene_Symbol = data_gene.name 
            if data_gene.mean( ) > 1 :
                is_normalized, scaleing_factor = True, data_gene.mean( )
                data_gene /= data_gene.mean( ) # normalize external data with mean value
            elif type( data_gene.values[ 0 ] ) is np.bool_ or type( data_gene.values[ 0 ] ) is  bool :
                data_gene = ( data_gene * series_organelle_amount - 0.05 ) * 1.05# if data_type of external data is boolean, hide data points where data = 0 by subtracting 0.1
        else : 
            is_proteome_data = True
            Gene_ID, Gene_Symbol = GET_Gene_ID_and_Symbol_from_input_gene( Gene )
            data_gene = df_all_input_genes.loc[ Gene_ID ]
        tumor_data_gene = data_gene[ list_tumor_samples_sorted ].values
        normal_data_gene = data_gene[ list_normal_samples_sorted ].values
        ax1.plot( tumor_pos, tumor_data_gene, 'o', color = T_gene_color, label = '(Tumor) ' + Gene_Symbol + is_proteome_data * ' Amount' )
        ax1.plot( normal_pos, normal_data_gene, 'o', color = N_gene_color, label = '(Normal) ' + Gene_Symbol + is_proteome_data * ' Amount' )
        if show_ref : # show scatter plots of data a given gene of reference samples if 'show_ref' is set to True
            ref_data_gene = data_gene[ list_ref_samples_sorted ].values
            ax1.plot( ref_pos, ref_data_gene, 'o', color = ref_gene_color, label = '(Ref) ' + Gene_Symbol + is_proteome_data * ' Amount' )
    organelle_min_amount, organelle_max_amount = series_organelle_amount.values.min( ), series_organelle_amount.values.max( )
    plt.xticks( np.arange( 0, organelle_max_amount, 0.1 ), np.round( np.arange( 0, organelle_max_amount, 0.1 ), 2 ) ) # set x_ticks 
    if is_normalized :
        plt.ylabel( Gene_Symbol ) # set label of gene_data as a y_label
        plt.yticks( np.arange( 0, 3, 0.2 ), np.round( np.arange( 0, 3, 0.2 ) * scaleing_factor, 2 ) ) # set y_ticks if gene_data has been normalized (priority to gene_data) 
    else :
        plt.ylabel( 'Distributions of ' + Gene_Set_description + ' Proteins' + plot_normalized * ' (Relative Amount to Pseudo Ref.)' )
    plt.xlim( left = series_organelle_amount.min( ) - 0.025, right = series_organelle_amount.max( ) + 0.025 ) # set x_range automatically
    plt.ylim( bottom = 0 ) # set y_limit
    plt.xlabel( Gene_Set_description + ' Amount From Relative Log Expression (RLE) Methods' + plot_normalized * ' Before Normalization' ) # annotate the plot
    plt.title( '{dataset_name} Proteomic data\n{description} Proteins Distribution (Total {n_genes} Proteins) {normalized}\n(Red = Tumor, Blue = Normal)'.format( dataset_name = dataset_name, description = Gene_Set_description, n_genes = str( num_Organelle_Genes ), normalized = plot_normalized * ' After RLE Normalization' ) )
    plt.legend( )
    plt.grid( False )


# In[ ]:


def PLOT_paired_T_vs_N( series_data, list_patients = None, alpha = 0.3, linewidth = 5, color_by_change = True ) :
    ''' Plot data in 'series_data' (where index is Sample_IDs of paired samples) on categorical plots (Normal, Tumor) so that change from
    tumor to normal can be tracked for individual patients with paired samples. Return fig and ax returned by plt.subplot '''
    if list_patients is None : # set default list_patients as all patients with paired samples
        list_patients = df_paired_samples.index.values
    else :
        list_patients = set( df_paired_samples.index ).intersection( list_patients )
    arr_paired_samples = df_paired_samples.loc[ list_patients ].values # retrive sample_IDs of each patient
    categories = [ 'Normal', 'Tumor' ] # define categories for plotting 
    fig, ax = plt.subplots( ) 
    for paired_samples_of_a_patient in arr_paired_samples : # plot tumor and normal data for each patient
        data_patient = series_data[ paired_samples_of_a_patient ]
        if color_by_change : # if 'color_by_change' is set to True, set color of plot according to Log2FC of Tumor to Normal
            color = cm.coolwarm( int( np.log2( data_patient[ 1 ] / data_patient[ 0 ] ) * 128 + 128 ) )
        else :
            color = 'k'
        ax.plot( categories, data_patient, '.-', lw = linewidth, color = color, alpha = alpha ) # overload alpha in color with 'alpha' argument
    ax.set_ylim( bottom = 0 ) # set y_limit
    return fig, ax # Return fig and ax returned by plt.subplot


# ### Basic functions for correlation analysis

# In[ ]:


def CORREL_s_df( s, df, print_message = True, correl_method = 'pearsonr', return_n_sample_with_zero_values = True ) :
    ''' Calculate Spearman correlation coefficient between series and df after alignment  '''
    l_indices_shared = list( set( df_count_norm.columns ).intersection( s.index.values ) )
    s = s.loc[ l_indices_shared ]
    df = df[ l_indices_shared ]
    data_s = s.values
    if correl_method == 'pearsonr' :
        df_result = pd.DataFrame( list( list( stats.pearsonr( data_s, data ) ) for data in df.values ), index = df.index.values, columns = [ 'correl_coeffi', 'p_value' ] )
    elif correl_method == 'kendalltau' :
        df_result = pd.DataFrame( list( list( stats.kendalltau( data_s, data ) ) for data in df.values ), index = df.index.values, columns = [ 'correl_coeffi', 'p_value' ] )
    else :
        df_result = pd.DataFrame( list( list( stats.spearmanr( data_s, data ) ) for data in df.values ), index = df.index.values, columns = [ 'correl_coeffi', 'p_value' ] )
    if return_n_sample_with_zero_values :
        df_result[ 'n_samples_with_zero_values' ] = ( df == 0 ).sum( axis = 1 )
    return df_result


# In[ ]:


def CORREL_gene_df( gene, df, df_target = None, print_message = True ) :
    ''' Calculate Spearman correlation coefficient (using NUMPY_spearman_correl) between given genes and other genes in the given df, and return the result as a pandas series
     if 'df_target' is given, perform correlation with genes in df_target  '''
    geneid, genesymbol = ( gene, gene ) if gene in df.index.values else GET_Gene_ID_and_Symbol_from_input_gene( gene ) # if a given label as a gene exist in df, assume df do not cotain genes.
    if geneid not in df.index.values :
        print( 'Gene does not exist in df' )
        return -1 
    s_gene = df.loc[ geneid ].dropna( )
    if print_message :
        print( 'Number of valid samples :', len( s_gene ) )
    df_target = df if df_target is None else df_target     
    s_gene, df_target_T = s_gene.align( df_target.T, 'inner' ) # retrive data of valid samples
    df_target = df_target_T.T.dropna( )  # drop NaN values before correlation
    data, geneids, data_gene = df_target.values, df_target.index.values, s_gene.values
    s = pd.Series( dict( ( geneid, NUMPY_spearman_correl( data_a_gene, data_gene ) ) for data_a_gene, geneid in zip( data, geneids ) ) )
    s.name = 'Correl_with_{}'.format( genesymbol )
    return s        


# In[ ]:


def CORREL_INTERNAL_Retrive_correl_coeff_from_scipy_modules( correl_result ) :
    '''  return correlation coefficient from the object scipy.stats functions and other functions are returning  '''
    if isinstance( correl_result, ( np.float64, float ) ) :
        return correl_result
    elif isinstance( correl_result, ( stats.stats.SpearmanrResult, stats.stats.KendalltauResult ) ) :
        return correl_result.correlation
    elif isinstance( correl_result, ( tuple ) ) :
        return correl_result[ 0 ]
    else :
        print( "invalid correlation result" )
        return -1


# In[ ]:


def CORREL_df( df, df_target = None, print_message = True, thres_minimum_n_samples = 10, axis = 1, method = 'kendalltau' ) :
    ''' Calculate Spearman correlation coefficient (using NUMPY_spearman_correl) between var of in the given df, and return the result as a pandas DataFrame (square matrix)
    df, df_target will be aligned on the given 'axis' (by dafault, aligned ) 
    if 'df_target' is given, calculate Spearman correlation between df and df_target (index = indices of df, columns = indices of df_target). DataFrame with NaN values are acceptable, but when
    number of valid samples for a var pair is smaller than 'thres_minimum_n_samples', data_value will be np.nan
    when there is no np.nan values in a dataframe and thus there is no need for masking np.nan values, scipy.stats methods are used for calculaating correlation coefficient and p_values
    'method' : a function name in scipy.stats module. if None is given, a default built-in custum function 'NUMPY_spearman_correl' will be used, which is significantly faster if number of observation is below 100 due to smaller overhead  '''
    flag_symmetric = True if df_target is None else False # set a flag to compute only half of a matrix if only one dataframe was given and thus correlation matrix is symmetric
    df_target = df if df_target is None else df_target # set df as df_target if only one dataframe was given 
    if axis == 0 : # if columns of the two dataframes become rows and columns of correlation matrix, transpose the two matrix 
        df, df_target = df.T, df_target.T
    df, df_target = df.align( df_target, join = 'inner', axis = 1 ) # align two dataframe before correlation analysis
    data_1, data_2 = df.values, df_target.values # retrive data of Normal ans Tumor samples
    mask_1, mask_2 = np.isnan( data_1 ), np.isnan( data_2 ) # get mask for NaN values
    l_var_1, n_var_1, l_var_2, n_var_2 = df.index.values, len( df ), df_target.index.values, len( df_target )
    arr_corr = np.zeros( ( n_var_1, n_var_2 ) ) # create 11 empty ( var x var ) numpy arrays filled zeros that will store computation results that can measure deviation from normal correlation
    Function_correl = stats.__dict__.get( method, NUMPY_spearman_correl ) # default correlation method is the default built-in custum function 'NUMPY_spearman_correl'
    if mask_1.sum( ) > 0 or mask_2.sum( ) > 0 :
        print( 'NaN values detected in aligned dataframes. Use mask to calculate correlation coefficient' )
        for index_var_1, data_var_1, mask_var_1 in zip( np.arange( len( data_1 ) ), data_1, mask_1 ) :
            for index_var_2, data_var_2, mask_var_2 in zip( np.arange( len( data_2 ) ), data_2, mask_2 ) :
                if flag_symmetric and index_var_2 > index_var_1 : 
                    continue
                mask_var_1_2 = ~ ( mask_var_1 | mask_var_2 ) # retrive mask for valid samples
                arr_corr[ index_var_1, index_var_2 ] = CORREL_INTERNAL_Retrive_correl_coeff_from_scipy_modules( Function_correl( data_var_1[ mask_var_1_2 ], data_var_2[ mask_var_1_2 ] ) ) if mask_var_1_2.sum( ) > thres_minimum_n_samples else np.nan
    else :
        for index_var_1, data_var_1 in zip( np.arange( len( data_1 ) ), data_1 ) :
            for index_var_2, data_var_2 in zip( np.arange( len( data_2 ) ), data_2 ) :
                if flag_symmetric and index_var_2 > index_var_1 :
                    continue
                else :
                    arr_corr[ index_var_1, index_var_2 ] = CORREL_INTERNAL_Retrive_correl_coeff_from_scipy_modules( Function_correl( data_var_1, data_var_2 ) ) # calculate a correlation coefficient                     
    df_corr = pd.DataFrame( arr_corr, index = l_var_1, columns = l_var_2 )
    df_corr.index.name, df_corr.columns.name = 'Indices of df', 'Indices of df_target'
    return PANDAS_relation_dataframe_make_symmetric( df_corr, symmetry_relationship = '=', diagonal_data = 1 ) if flag_symmetric else df_corr


# In[ ]:


def CORREL_s_with_df( s, df, print_message = True ) :
    ''' Calculate Spearman correlation coefficient (using NUMPY_spearman_correl) between given genes and other genes in the given df, and return the result as a pandas series
     if 'df_target' is given, perform correlation with genes in df_target  '''
    s, dfT = s.align( df.T, 'inner' ) # retrive data of valid samples
    df = dfT.T.dropna( ) # drop NaN values before correlation
    if print_message : print( "number of features used for correlation analysis: {}".format( len( df.columns ) ) )
    data_df, arr_index_df, data_s = df.values, df.index.values, s.values
    s_correl_result = pd.Series( dict( ( index, NUMPY_spearman_correl( data_of_an_index, data_s ) ) for data_of_an_index, index in zip( data_df, arr_index_df ) ) )
    s_correl_result.name = 'CORREL_s_with_df'
    return s_correl_result


# ### Functions for GSEA

# In[ ]:


def GSEA_RESULT_Search_index( df, query ) :
    return df[ list( True if query in index else False for index in df.index.values ) ]


# ### Functions During Analysis

# In[ ]:


"""def CACULATE_correlation__sample_set__Gene( sample_set, gene = None, correl_type = 'pearson', drop_na = True ) :
    '''
    perform correlation with a given gene in a given set of samples
    return DataFrame of correlation result with a given gene 
    '''
    if drop_na : # drop NaN values in df_proteome to increase the accuracy of correlation analysis
        df_selected = df_proteome[ sample_set ].dropna( ) # select only given samples in proteome dataset without NaN values
    else :
        df_selected = df_proteome[ sample_set ] # select only given samples in proteome dataset
    gene_id = Gene_2_Gene_ID( gene ) # retrive gene_id
    if gene_id == -1 : # if gene is invalid, end the method
        return
    target_data = df_selected.loc[ gene_id ] # retrive data of target gene, or "bait" gene 
    correlation = dict( correl_coef = [ ], Gene_Symbol = [ ], p_value = [ ], Gene_ID = [ ], Gene_Name = [ ] )
    if correl_type == 'pearson' : # perform correlation analysis
        for Gene_ID in df_selected.index :
            data = df_selected.loc[ Gene_ID ]
            result = stats.pearsonr( data, target_data )
            correlation[ 'correl_coef' ].append( result[ 0 ] )
            correlation[ 'p_value' ].append( result[ 1 ] )
            correlation[ 'Gene_ID' ].append( Gene_ID ) 
            correlation[ 'Gene_Name' ].append( dict_ID_Symbol[ Gene_ID ][ 1 ] ) 
            correlation[ 'Gene_Symbol' ].append( dict_ID_Symbol_simple[ Gene_ID ] ) 
    else : # if correlation type is spearman's correlation
        for Gene_ID in df_selected.index :
            data = df_selected.loc[ Gene_ID ]
            result = stats.spearmanr( data, target_data, nan_policy = 'omit' )
            correlation[ 'correl_coef' ].append( result.correlation )
            correlation[ 'p_value' ].append( result.pvalue )
            correlation[ 'Gene_ID' ].append( Gene_ID ) 
            correlation[ 'Gene_Name' ].append( dict_ID_Symbol[ Gene_ID ][ 1 ] ) 
            correlation[ 'Gene_Symbol' ].append( dict_ID_Symbol_simple[ Gene_ID ] ) 
    # convert correlation result into DataFrame and handle DataFrame
    df = pd.DataFrame( correlation )
    df = df.set_index( 'Gene_ID' )
    df = df.dropna( )
    df = df.sort_values( [ 'correl_coef' ], axis = 'index' )
    return df"""


# In[ ]:


def CACULATE_correlation__sample_list__Gene( df, sample_list, gene = None, correl_type = 'pearson', drop_na = True ) :
    '''    perform correlation with a given gene in a given set of samples.    return DataFrame of correlation result with a given gene.    '''
    if drop_na : # drop NaN values in df_proteome to increase the accuracy of correlation analysis
        df_selected = df_proteome[ sample_list ].dropna( ) # select only given samples in proteome dataset without NaN values
    else :
        df_selected = df_proteome[ sample_list ] # select only given samples in proteome dataset
    gene_id = Gene_2_Gene_ID( gene ) # retrive gene_id
    if gene_id == -1 : # if gene is invalid, end the method
        return
    target_data = df_selected.loc[ gene_id ] # retrive data of target gene, or "bait" gene 
    correlation = dict( correl_coef = [ ], Gene_Symbol = [ ], p_value = [ ], Gene_ID = [ ], Gene_Name = [ ] )
    if correl_type == 'pearson' : # perform correlation analysis
        for Gene_ID in df_selected.index :
            data = df_selected.loc[ Gene_ID ]
            result = stats.pearsonr( data, target_data )
            correlation[ 'correl_coef' ].append( result[ 0 ] )
            correlation[ 'p_value' ].append( result[ 1 ] )
            correlation[ 'Gene_ID' ].append( Gene_ID ) 
            correlation[ 'Gene_Name' ].append( dict_ID_Symbol[ Gene_ID ][ 1 ] ) 
            correlation[ 'Gene_Symbol' ].append( dict_ID_Symbol_simple[ Gene_ID ] ) 
    else : # if correlation type is spearman's correlation
        for Gene_ID in df_selected.index :
            data = df_selected.loc[ Gene_ID ]
            result = stats.spearmanr( data, target_data, nan_policy = 'omit' )
            correlation[ 'correl_coef' ].append( result.correlation )
            correlation[ 'p_value' ].append( result.pvalue )
            correlation[ 'Gene_ID' ].append( Gene_ID ) 
            correlation[ 'Gene_Name' ].append( dict_ID_Symbol[ Gene_ID ][ 1 ] ) 
            correlation[ 'Gene_Symbol' ].append( dict_ID_Symbol_simple[ Gene_ID ] ) 
    # convert correlation result into DataFrame and handle DataFrame
    df = pd.DataFrame( correlation )
    df = df.set_index( 'Gene_ID' )
    df = df.dropna( )
    df = df.sort_values( [ 'correl_coef' ], axis = 'index' )
    return df


# In[ ]:


def Shrink_by_p_value( p_val, p_value_filter = 0.25 ):
    '''    Calculate weight of an entry based on the p-value and p_value_filter.    p_val should be lower than 1 and larger than 0    '''
    if np.isnan( p_val ) :
        return np.nan
    ratio = np.log( p_value_filter * 1.5 ) / np.log( p_val ) # calculate the ratio of the log of p-value filter to p-value
    return 0 if ratio > 1 else ( 1 - ratio ) # return the weight of this entry calcluated from p-value


# In[ ]:


"""def p_value_scoring_ratio_large_slope( p_val, p_value_filter = 0.25 ) : # a function for Log2FC shrinkage
    '''    Calculate weight of an entry based on the p-value and p_value_filter
    p_val should be lower than 1 and larger than 0    '''
    if np.isnan( p_val ) :
        return np.nan
    ratio = np.log( p_value_filter * 1.5 ) / np.log( p_val ) # calculate the ratio of the log of p-value filter to p-value
    return 0 if ratio > 1 else ( 1 - ratio ) # return the weight of this entry calcluated from p-value"""


# In[ ]:


def Export_Log2FC_tumor_vs_normal_Pathview( df, name_of_df = '' ) :
    '''
    Export Log2FC tumor_vs_normal csv file for Pathview visualization in the result_folder
    Log2FC is adjusted by corresponding t-test p_value to decrease rate of false changes
    '''
    Gene_IDs = df.index.values # retrive Gene_IDs 
    normal_data = df[ dict_type_sample[ 'normal' ] ].values # retrive tumor and normal data 
    tumor_data = df[ dict_type_sample[ 'tumor' ] ].values
    normal_data = np.ma.masked_array( normal_data, np.isnan( normal_data ) ) # apply mask on np.nan values of tumor and normal data 
    tumor_data = np.ma.masked_array( tumor_data, np.isnan( tumor_data ) )

    p_values = [ ] # calculate p_values and adjust Log2FC values according to the p_value
    Log2FC_tumor_by_normal = [ ]
    for intIndex in range( len( Gene_IDs ) ) : # for each entry
        normal = normal_data[ intIndex ] # retrive normal and tumor data for the entry
        tumor = tumor_data[ intIndex ]
        p_value = stats.ttest_ind( normal, tumor, nan_policy = 'omit' ).pvalue # perform t-test (with ignoring np.nan values) and add p_value to the list
        p_values.append( p_value ) 
        # calculate log2 fold change and multiply weight based on p_value, and add Log2FC to the list
        Log2FC_tumor_by_normal.append( np.log2( np.average( tumor ) / np.average( normal ) ) * p_value_scoring_ratio_large_slope( p_value, p_value_filter = 0.001 ) )

    # create dataframe from t-test, log2FC values
    dict_normal_vs_tumor_FC_proteome = dict( Gene_ID = Gene_IDs, Log2FC_tumor_by_normal = Log2FC_tumor_by_normal, p_value = p_values )
    df = pd.DataFrame( dict_normal_vs_tumor_FC_proteome )
    df = df.set_index( 'Gene_ID' )
    df = df.sort_values( [ 'Log2FC_tumor_by_normal' ], axis = 'index' ) # sort the result for visualization
    #invalid_Indices = df[ pd.isnull( df.p_value ) ].index.values # retrive invalid entries by identifing invalid p_values
    #df = df.drop( labels = invalid_Indices, axis = 'index' ) # drip entries with invalid values
    df_to_export = df.drop( labels = 'p_value', axis = 'columns' ) # drop p_value colume for the use in Pathview R package
    df_to_export.to_csv( result_folder + 'normal_vs_tumor_FC_proteome' + name_of_df + '.csv' ) # export dataframe to csv


# In[ ]:


def normalization( df, gene_set ) :
    ''' normalize protein amounts in df by using average protein amounts of a given gene_set, and return it as a dataframe '''
    gene_set = set( gene_set ).intersection( set( df.index.values ) )
    average_gene_set = np.average( df.loc[ gene_set ].values, axis = 0 )
    df_normalization = deepcopy( df )
    df_normalization.iloc[ :, : ] = df_normalization.values / average_gene_set
    return df_normalization


# In[ ]:


def GET_dataframe_lower_triangle( df ) :
    ''' perform np.tril on DataFrame values and return a copy of DataFrame that only containing lower triangle (np.tril) '''
    df = deepcopy( df )
    df.iloc[ :, : ] = np.tril( df )
    return df


# In[ ]:


def LINREGRESS_make_dataframe_symmetric( df, mask_xy_switched_tril, mask_tril, fill_diagnal_value, method_xy_switch = lambda a, b : a ) :
    ''' a fuction intended to be used in linear regression analysis for selecting values using given masks and makes symmetrix matrix, while filling diagnal entries with 'fill_diagnal_value' 
    'method_xy_switch' is used for calculating values when xy axes are switched (for example, df_slope needs lambda a, b : 1 / a ) 'a' = df, 'b' = a flag that is True when 'a' is df_xy_switched and not selected, False when 'a' is not switched and selected '''
    df_xy_switched_tril, df_tril = GET_dataframe_lower_triangle( method_xy_switch( df.transpose( ), True ) ), GET_dataframe_lower_triangle( df ) # retrive inverse (1/s) of slopes where x and y are switched  # remove data above the diagonal line to simplify conputation
    df_selected_tril = df_xy_switched_tril * mask_xy_switched_tril + df_tril * mask_tril    
    df_selected_triu = deepcopy( df_selected_tril )
    df_selected_triu.values[ indices_triupper ] = 1
    df_selected_triu = method_xy_switch( df_selected_triu.transpose( ), False )
    df_selected_triu.values[ indices_trilower ] = 0
    df_selected = df_selected_tril + df_selected_triu
    df_selected.values[ indices_diagonal ] = fill_diagnal_value
    return df_selected


# In[ ]:


def LINREGRESS_intercept_xy_switch_method( df, is_xy_switched_and_not_selected ) :
    if is_xy_switched_and_not_selected :
        return - ( df.transpose( ) / df_slope ).transpose( )
    else :
        return - ( df.transpose( ) / df_slope_selected ).transpose( )


# In[ ]:


def mstats_linregress_theilsopes_average( a, b ) :
    ''' calculate average slope and intercept of a (x-axis), b (y-axis) using both linregress and theilsopes functions in mstats modules (can handle masked arrays) '''
    slope, intercept, rvalue, pvalue, stderr = mstats.linregress( a, b )
    s_t, i_t, _, _ = mstats.theilslopes( b, a )
    slope, intercept = ( slope + s_t ) / 2, ( intercept + i_t ) / 2
    return slope, intercept, rvalue, pvalue, stderr
def stats_linregress_theilsopes_average( a, b ) :
    ''' calculate average slope and intercept of a (x-axis), b (y-axis) using both linregress and theilsopes functions in stats modules  '''
    slope, intercept, rvalue, pvalue, stderr = stats.linregress( a, b )
    s_t, i_t, _, _ = stats.theilslopes( b, a )
    slope, intercept = ( slope + s_t ) / 2, ( intercept + i_t ) / 2
    return slope, intercept, rvalue, pvalue, stderr


# ### Functions for finding an elbow point

# In[ ]:


def ELBOW_Find( l, sort_array = True, return_index_of_an_elbow_point = False ) :
    """ # 2020-12-22 23:45:13 
    Find an elbow point based on the rotation method.
    (return the value at the elbow point)
    'return_index_of_an_elbow_point' : if True, return an index of an elbow point in a sorted array (ascending = True). if False, return the value of an elbow point 
    """
    # preprocess input array
    if isinstance( l, ( pd.Series ) ) :
        l = l.values
    if not isinstance( l, ( np.ndarray ) ) :
        l = np.array( l, dtype = float ) # convert data to array with floating point numbers
    LIST_COUNT
    if sort_array :
        l = np.sort( l ) # sort the given data
    width_y = l[ -1 ] - l[ 0 ] # get the difference between max and min values of the given data values (y)
    x = np.arange( len( l ) ) / len( l ) * width_y # get scaled x-coordinates of the data (scaled floating point indices of the sorted data array) (scaled to the difference between min and max values in the data array)
#     sine_theta = width_y / ( ( l[ - 1 ] - l[ 0 ] ) ** 2 + ( x[ - 1 ] - x[ 0 ] ) ** 2 ) ** 0.5
    sin_45_degree = 1 / 2 ** 0.5
    arr_rotation = np.array( [ [ sin_45_degree, sin_45_degree ], [ - sin_45_degree, sin_45_degree ] ] )
    x_rotated, y_rotated = np.matmul( arr_rotation, np.vstack( [ x, l ] ) )
    index_elbow_point = y_rotated.argmin( )
    
    return index_elbow_point if return_index_of_an_elbow_point else l[ index_elbow_point ] # return the index or the value of an elbow point


# ### Functions using Multiprocessing librarys

# In[ ]:


def Multiprocessing( arr, Function, n_threads = 12, dir_temp = '/tmp/', Function_PostProcessing = None ) : 
    """ # 2020-12-15 22:41:06 
    split inputs given by 'arr' into 'n_threads' number of temporary files, and folks 'n_threads' number of processes running a function given by 'Function' by givning a directory of each temporary file as an argument. if arr is DataFrame, the temporary file will be split DataFrame (tsv format) with column names, and if arr is 1d or 2d array, the temporary file will be tsv file without header 
    
    'Function_PostProcessing' : if given, Run the function before removing temporary files at the given temp folder. uuid of the current session and directory of the temporary folder are given as arguments to the function.
    """
    if isinstance( arr, ( list ) ) : # if a list is given, convert the list into a numpy array
        arr = np.array( arr )
    str_uuid = UUID( ) # create a identifier for making temporary files
    l_dir_file = list( ) # split inputs given by 'arr' into 'n_threads' number of temporary files
    if isinstance( arr, pd.DataFrame ) : # if arr is DataFrame, the temporary file will be split DataFrame (tsv format) with column names
        for index_chunk in range( n_threads ) :
            dir_file_temp = dir_temp + str_uuid + '_' + str( index_chunk ) + '.tsv.gz'
            arr.iloc[ index_chunk : : n_threads ].to_csv( dir_file_temp, sep = '\t', index = False )
            l_dir_file.append( dir_file_temp )
    else : # if arr is 1d or 2d array, the temporary file will be tsv file without header
        l_chunk = LIST_Split( arr, n_threads )
        for index, arr in enumerate( l_chunk ) : # save temporary files containing inputs
            dir_file_temp = dir_temp + str_uuid + '_' + str( index ) + '.tsv'
            if len( arr.shape ) == 1 : df = pd.DataFrame( arr.reshape( arr.shape[ 0 ], 1 ) )
            elif len( arr.shape ) == 2 : df = pd.DataFrame( arr )
            else : print( 'invalid inputs: input array should be 1D or 2D' ); return -1
            df.to_csv( dir_file_temp, sep = '\t', header = None, index = False )
            l_dir_file.append( dir_file_temp )

    with Pool( n_threads ) as p : l = p.map( Function, l_dir_file ) # use multiple process to run the given function
        
    if Function_PostProcessing is not None :
        Function_PostProcessing( str_uuid, dir_temp ) 
        
    for dir_file in glob.glob( dir_temp + str_uuid + '*' ) : os.remove( dir_file ) # remove temporary files
    return l # return mapped results


# ### Numpy Utility Functions

# In[ ]:


def NUMPY_Calculate_Average( arr ) : # 2020-07-30 16:17:33 
    ''' For a given np.array with float dtype, count number of non-nan values, and return the average of non-nan values along with the count of non-nan values. if there is no non-nan values, return np.nan and 0 ''' 
    arr_valid = arr[ ~ np.isnan( arr ) ]
    n_non_nan_values = len( arr_valid )
    avg = arr_valid.sum( ) / n_non_nan_values if n_non_nan_values > 10 else ( sum( arr_valid ) / n_non_nan_values if n_non_nan_values > 0 else np.nan )
    return avg, n_non_nan_values


# In[ ]:


def LISTNUMPY_List_comprehension_using_numpy( arr, f_map_return, dtype, f_map_if = None, f_map_else_return = None, f_if = None ) :
    ''' Numpy implementation of list comprehension, which takes less memory than list comprehension, which treat every element as an object and thus takes more space '''
    if f_if is not None : #  if 'f_if' is not None, filter out entries that makes 'f_if' return False
        arr_mask = np.zeros_like( arr, dtype = bool )
        for index in np.arange( len( arr ) ) :
            arr_mask[ index ] = f_if( arr[ index ] )
        arr = arr[ arr_mask ]
    arr_mapped = np.zeros_like( arr, dtype = dtype )
    if f_map_if is None :
        for index in np.arange( len( arr ) ) :
            arr_mapped[ index ] = f_map_return( arr[ index ] )
    else :
        if f_map_else_return is None : # if f_map_else_return is None, return np.nan if 'f_map_if' returns False for the entry
            for index in np.arange( len( arr ) ) :
                arr_mapped[ index ] = f_map_return( arr[ index ] ) if f_map_if( arr[ index ] ) else np.nan
        else :
            for index in np.arange( len( arr ) ) :
                arr_mapped[ index ] = f_map_return( arr[ index ] ) if f_map_if( arr[ index ] ) else f_map_else_return( arr[ index ] )
    return arr_mapped


# In[ ]:


def NUMPY_POLY_Get_y( p, x ) :
    ''' p : p(x) = p[0] * x**deg + ... + p[deg], infer degree of polynomial from length of p  '''
    deg = len( p ) - 1
    y = 0.0
    for index in np.arange( deg + 1 ) :
        y += p[ index ] * x ** ( deg - index )
    return y


# In[ ]:


def NUMPY_UTIL_Does_contain_NaN( arr, print_message = True ) :
    ''' Return True if a given array contains np.nan, regardless of dtype '''
    n_NaN = 0
    for entry in arr :
        if isinstance( entry, ( float, np.float64, np.float32 ) ) and np.isnan( entry ) :
            n_NaN += 1
    if n_NaN > 0 :
        if print_message :
            print( "[NUMPY_UTIL_Does_contain_NaN] The given array contains {} n_NaN values".format( n_NaN ) )
        return True
    else :
        return False


# ### Optimized Functions (custom implementation of Scipy and numpy methods)

# In[ ]:


def NUMPY_Log_Average( arr ):
    return np.exp( np.log( arr ).mean( ) )


# In[ ]:


def NUMPY_GET_mean_std( data ) :
    ''' calculate mean and standard deviation of a given array of data using np.dot, which involves vector calculation (very fast). when number of elements is small (~100), this method is 
    3 times faster than applying data.mean() and data.std() separately.  calculate sample standard deviation (divided by N - 1 instead of N). return  mean, std     '''
    N = data.size
    mean = np.dot( data, np.ones( N ) ) / N
    deviation = data - mean
    std = np.sqrt( np.dot( deviation, deviation ) / ( N - 1 ) ) # calculate sample standard deviation (degree of freedom = 1)
    return mean, std


# In[ ]:


def NUMPY_spearman_correl( data_1, data_2 ) :
    ''' Use np.argsort() and np.dot() to reduce computation time. When number of elements is small (~100), it is 50 times faster than scipy.stats.spearmanr method 
    Since np.argsort() give identical values distinct rank, this method can give incorrect spearman correlation coefficient if given data have identical values (for example, multiple zeros)  '''
    N = data_1.size
    rank_1 = data_1.argsort().argsort() # retrive a rank of given data, assuming all there are no duplicates in the data
    rank_2 = data_2.argsort().argsort()
    diff_of_rank = rank_1 - rank_2
    return 1 - 6 * np.dot( diff_of_rank, diff_of_rank ) / ( N * ( N ** 2 - 1 ) ) # calculate spearman's rank correlation 


# In[ ]:


def NUMPY_least_square_linear_regress( x, y ) :
    ''' Perform simple linear regression (least square linear regression), which is implemented in scipy.stats.linregress. When number of elements is moderate ( < 10,000 ), 
    this method is 10 ~ 20 times faster than linregress method. (when n_elements = 100, its 20 times faster than linregress) '''
    N = x.size
    arr_ones = np.ones( N )
    sum_x, sum_y, sum_x2, sum_xy = np.dot( arr_ones, x ), np.dot( arr_ones, y ), np.dot( x, x ), np.dot( x, y ) # calculate sum of x, y, x^2, xy that will be used to calculate the slope 
    slope = ( N * sum_xy - sum_x * sum_y ) / ( N * sum_x2 - sum_x ** 2 )
    return slope, ( sum_y - slope * sum_x ) / N # return slope and intercept (intercept calculated using mean of x and y)


# In[ ]:


def skip_diag_masking( A ) :
    ''' remove diagnoal entries and reshape the matrix (example, 3 x 3 matrix becomes 3 x 2 matrix) '''
    return A[ ~ np.eye( A.shape[ 0 ],dtype = bool ) ].reshape( A.shape[ 0 ], - 1 )


# In[ ]:


def NUMPY_theil_sen_regression( data_1, data_2 ) :
    ''' brute force (instead of random sampling) calculation of Theil-Sen linear regression slope. For small number of elements (<100), its about 2 times faster than its Scipy implementation '''
    N = data_1.size
    arr_x, arr_y = np.zeros( ( 2, N, N ) ) # create arrays to which x and y coordinates are broadcasted (and duplicated N times) 
    arr_x[ :, : ], arr_y[ :, : ] = data_1, data_2 # broad cast x and y data to N x N matrix
    slope = np.median( skip_diag_masking( arr_y.T - data_2 ) / skip_diag_masking( arr_x.T - data_1 ) ) # calculate median of slopes of all pairs of distinct points (diagonals are removed since they are zeros (slope between identical points) and do not have slopes)
    return slope, np.median( data_2 ) - slope * np.median( data_1 ) # return slope and an intercept


# In[ ]:


def NUMPY_LIST_argmax( a_list ):
    ''' Receive a list and return index of the maximum value. Due to hugh overhead of np.argmax, this function is 12 times faster when there are only two elements '''
    return a_list.index( max( a_list ) )


# In[ ]:


def NUMPY_LIST_argmin( a_list ):
    ''' Receive a list and return index of the minimum value. Due to hugh overhead of np.argmin, this function is 12 times faster when there are only two elements '''
    return a_list.index( min( a_list ) )


# ##### Functions not for optimization but primarily use numpy functions

# In[ ]:


def NUMPY_expand_narrow_arr( df, width_column = None ) :
    if width_column is None :
        width_column = int( len( df ) / 10 )
    elif width_column == 'fit' :
        width_column = int( len( df ) / len( df.T ) )
    data = df.T.values
    arr_expanded = np.zeros( ( 0, len( df ) ) )
    for data_col in data :
        col_expanded = np.zeros( ( width_column, len( df ) ) )
        col_expanded[ :, : ] = data_col
        arr_expanded = np.vstack( ( arr_expanded, col_expanded ) )
    return arr_expanded.T


# ### Functions for Linear Regression Analysis (DEVIATION_FROM_NORMAL_CORRELATION Analysis)

# In[ ]:


arr_halfs = np.ones( 2 ) / 2 # define an object that is required for a method that uses below method 
def LINREGRESS_OPTIMIZE_calculate_std_of_rejections( slope, x_y_centered ) :
    ''' A function that rotate centered data according from a given slope to x_axis, and calculate the standard deviation of rejections (y-axis after rotation)  '''
    return NUMPY_GET_mean_std( np.matmul( np.array( [ [ 1, slope ], [ - slope, 1 ] ] ) / np.sqrt( 1 + slope ** 2 ), x_y_centered )[ 1 ] )[ 1 ] # create a rotation matrix, perform rotation, select y values, and return std of the values


# In[ ]:


an_intercept_vector = np.array( [ [ 0.0 ], [ 0.0 ] ] ) # an empty 2D vector then will store the value of an intercept and become an intercept vector 
def LINREGRESS_Projection_Based_Least_Square_Linear_Regression( gene_1_data, gene_2_data, gene_1_2 ) :
    ''' Perform Projection-Based Least Square Linear Regression on a given set of data (gene_1_data = x_data, gene_2_data = y_data, gene_1_2 = stacked gene_1 and gene_2 data) 
    require an_intercept_vector = np.array( [ [ 0 ], [ 0 ] ] ) defined outside this function, which can speed up the calculation
    return optimal_slope and optimal_intercept '''
    n = gene_1_data.size
    S_x, S_y = gene_1_data.sum( ), gene_2_data.sum( )
    S_x_y = ( gene_1_data * gene_2_data ).sum( )
    S_x2, S_y2 = ( gene_1_data * gene_1_data ).sum( ), ( gene_2_data * gene_2_data ).sum( )
    minus_b_by_2a = - ( S_y ** 2 - S_x ** 2 + n * ( S_x2 - S_y2 ) ) / ( ( n * S_x_y - S_x * S_y ) * 2 )
    sqrt_value = np.sqrt( 1 + minus_b_by_2a ** 2 ) # inside square root, the value is always larger than 1 
    list_slopes = [ minus_b_by_2a - sqrt_value, minus_b_by_2a + sqrt_value ] # calculate the two possible optimal slopes and intercepts
    list_intercepts = [ ( S_y - S_x * list_slopes[ 0 ] ) / n, ( S_y - S_x * list_slopes[ 1 ] ) / n ]
    an_intercept_vector[ 1 ] = list_intercepts[ 0 ]
    std_rejection_1 = LINREGRESS_OPTIMIZE_calculate_std_of_rejections( list_slopes[ 0 ], gene_1_2 - an_intercept_vector )
    an_intercept_vector[ 1 ] = list_intercepts[ 1 ]
    std_rejection_2 = LINREGRESS_OPTIMIZE_calculate_std_of_rejections( list_slopes[ 1 ], gene_1_2 - an_intercept_vector )
    index_min = NUMPY_LIST_argmin( [ std_rejection_1, std_rejection_2 ] )
    return list_slopes[ index_min ], list_intercepts[ index_min ]


# In[ ]:


def DEVIATION_FROM_NORMAL_CORRELATION_calculation( df, Gene_Set = None, n_std_for_outliers = 3.5 ) :
    ''' For each gene pair in Gene_Set (default : all genes in df), perform projection-based least square linear regression for normal data, while masking outliers outside 'n_std_for_outliers'
    from the mean, and calculate mean rejection (Tumor data only since mean rejection in Normal data is 0), standard deviation(std) of rejection, std of projection, correlation between rejections 
    and projections, correlation before projection (Correlation Matrix of df with Spearman's coefficient), and lastly, deviation from origin from the slope and intercept (same sign with intercept). '''
    if Gene_Set is not None :
        df = PANDAS_Subset( df, Gene_Set )
    n_genes = len( df )
    data = df.values # row = gene, col = sample
    Gene_IDs = df.index.values

    T_data, N_data = df[ T_samples ].values, df[ N_samples ].values # retrive data of Normal ans Tumor samples
    T_mask, N_mask = OUTLIERS_GET_mask_for_outliers( T_data, n_std_for_outliers = n_std_for_outliers ), OUTLIERS_GET_mask_for_outliers( N_data, n_std_for_outliers = n_std_for_outliers ) # get mask for outliers
    arr_T_corr, arr_N_corr, arr_slopes, arr_intercepts, arr_T_rej_mean, arr_T_rej_std, arr_T_proj_std, arr_T_dev_corr, arr_N_rej_std, arr_N_proj_std, arr_N_dev_corr = np.zeros( ( 11, n_genes, n_genes ) ) # create 11 empty ( n_genes x n_genes ) numpy arrays filled zeros that will store computation results that can measure deviation from normal correlation

    for index_gene_1, T_data_gene_1, T_mask_gene_1, N_data_gene_1, N_mask_gene_1 in zip( np.arange( n_genes ), T_data, T_mask, N_data, N_mask ) :
        for index_gene_2, T_data_gene_2, T_mask_gene_2, N_data_gene_2, N_mask_gene_2 in zip( np.arange( n_genes ), T_data, T_mask, N_data, N_mask ) :
            if index_gene_2 > index_gene_1 :
                continue
            T_mask_gene_1_2 = ~ ( T_mask_gene_1 | T_mask_gene_2 ) # retrive a mask for non_outliers in the pair of datavalues of genes
            N_mask_gene_1_2 = ~ ( N_mask_gene_1 | N_mask_gene_2 )
            T_data_gene_1_no_outliers, T_data_gene_2_no_outliers = T_data_gene_1[ T_mask_gene_1_2 ], T_data_gene_2[ T_mask_gene_1_2 ] # retrive non_outliers data
            N_data_gene_1_no_outliers, N_data_gene_2_no_outliers = N_data_gene_1[ N_mask_gene_1_2 ], N_data_gene_2[ N_mask_gene_1_2 ]
            T_data_gene_1_2 = np.vstack( ( T_data_gene_1_no_outliers, T_data_gene_2_no_outliers ) )  # create an array of vectors of postitions (gene_1 = x axis, gene_2 = y axis)
            N_data_gene_1_2 = np.vstack( ( N_data_gene_1_no_outliers, N_data_gene_2_no_outliers ) )
            slope, intercept = LINREGRESS_Projection_Based_Least_Square_Linear_Regression( N_data_gene_1_no_outliers, N_data_gene_2_no_outliers, N_data_gene_1_2 ) # retrive slope and intercept of the gene_pair
            arr_slopes[ index_gene_1, index_gene_2 ], arr_intercepts[ index_gene_1, index_gene_2 ] = slope, intercept # put calculated index and slope into the array
            an_intercept_vector[ 1 ] = intercept # define vector of position of intercept
            R = np.array( [ [ 1, slope ], [ - slope, 1 ] ] ) / np.sqrt( 1 + slope ** 2 ) # calculation rotation matrix that rotate linear regression line to x-axis
            T_data_gene_1_2_projected = np.matmul( R, T_data_gene_1_2 - an_intercept_vector ) # subtract intercept_vector and perform projection to the linear regression line using the rotation matrix.
            N_data_gene_1_2_projected = np.matmul( R, N_data_gene_1_2 - an_intercept_vector )
            arr_T_rej_mean[ index_gene_1, index_gene_2 ], arr_T_rej_std[ index_gene_1, index_gene_2 ] = NUMPY_GET_mean_std( T_data_gene_1_2_projected[ 1 ] ) # Than, calculating mean and standard deviation of rejections of the points 
            _, arr_N_rej_std[ index_gene_1, index_gene_2 ] = NUMPY_GET_mean_std( N_data_gene_1_2_projected[ 1 ] )
            _, arr_T_proj_std[ index_gene_1, index_gene_2 ] = NUMPY_GET_mean_std( T_data_gene_1_2_projected[ 0 ] ) # Than, calculating mean and standard deviation of rejections of the points 
            _, arr_N_proj_std[ index_gene_1, index_gene_2 ] = NUMPY_GET_mean_std( N_data_gene_1_2_projected[ 0 ] )
            arr_T_dev_corr[ index_gene_1, index_gene_2 ] = NUMPY_spearman_correl( T_data_gene_1_2_projected[ 0 ], T_data_gene_1_2_projected[ 1 ] ) # calculate Spearman's rank correlation between projections and rejections of points
            arr_N_dev_corr[ index_gene_1, index_gene_2 ] = NUMPY_spearman_correl( N_data_gene_1_2_projected[ 0 ], N_data_gene_1_2_projected[ 1 ] )
            arr_T_corr[ index_gene_1, index_gene_2 ] = NUMPY_spearman_correl( T_data_gene_1_no_outliers, T_data_gene_2_no_outliers ) # calculate Spearman's rank correlation before projection
            arr_N_corr[ index_gene_1, index_gene_2 ] = NUMPY_spearman_correl( N_data_gene_1_no_outliers, N_data_gene_2_no_outliers )       

    df_T_corr, df_N_corr, df_slopes, df_intercepts, df_T_rej_mean, df_T_rej_std, df_T_proj_std, df_T_dev_corr, df_N_rej_std, df_N_proj_std, df_N_dev_corr = pd.DataFrame( arr_T_corr, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_N_corr, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_slopes, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_intercepts, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_T_rej_mean, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_T_rej_std, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_T_proj_std, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_T_dev_corr, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_N_rej_std, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_N_proj_std, index = Gene_IDs, columns = Gene_IDs ), pd.DataFrame( arr_N_dev_corr, index = Gene_IDs, columns = Gene_IDs ) # convert numpy arrays into DataFrames
    df_slopes = PANDAS_relation_dataframe_make_symmetric( df_slopes, symmetry_relationship = '1/' )
    df_intercepts = df_intercepts + ( - df_intercepts / df_slopes ).T # calculate half of intercept data from the other half with slope data
    df_T_rej_mean = PANDAS_relation_dataframe_make_symmetric( df_T_rej_mean, symmetry_relationship = '-' )
    df_T_rej_std, df_N_rej_std = PANDAS_relation_dataframe_make_symmetric( df_T_rej_std, symmetry_relationship = '=' ), PANDAS_relation_dataframe_make_symmetric( df_N_rej_std, symmetry_relationship = '=' )
    df_T_proj_std, df_N_proj_std = PANDAS_relation_dataframe_make_symmetric( df_T_proj_std, symmetry_relationship = '=' ), PANDAS_relation_dataframe_make_symmetric( df_N_proj_std, symmetry_relationship = '=' )
    df_T_dev_corr, df_N_dev_corr = PANDAS_relation_dataframe_make_symmetric( df_T_dev_corr, symmetry_relationship = '-', diagonal_data = 0 ), PANDAS_relation_dataframe_make_symmetric( df_N_dev_corr, symmetry_relationship = '-', diagonal_data = 0 )
    df_T_corr, df_N_corr = PANDAS_relation_dataframe_make_symmetric( df_T_corr, symmetry_relationship = '=', diagonal_data = 1 ), PANDAS_relation_dataframe_make_symmetric( df_N_corr, symmetry_relationship = '=', diagonal_data = 1 )
    df_dev_origin = df_intercepts / np.sqrt( df_slopes.values ** 2 + 1 ) # calculate deviation from origin (sign is same with intercept) from df_intercepts and df_slopes. deviation from origin = i / sqrt( s**2 + 1 )
    return dict( df_T_corr = df_T_corr, df_N_corr = df_N_corr, df_slopes = df_slopes, df_intercepts = df_intercepts, df_T_rej_mean = df_T_rej_mean, df_T_rej_std = df_T_rej_std, df_T_proj_std = df_T_proj_std, df_T_dev_corr = df_T_dev_corr, df_N_rej_std = df_N_rej_std, df_N_proj_std = df_N_proj_std, df_N_dev_corr = df_N_dev_corr, df_dev_origin = df_dev_origin )


# In[ ]:


def DEVIATION_NORMAL_CORR_a_gene_pair_with_switching( gene_1, gene_2, df = None, n_std_for_outliers = 3.5, mask_zero_with_nan = True ) :
    ''' calculates 14 metrics for projection-based linear regression analysis for a given pair of genes on 'df' data using 'DEVIATION_NORMAL_CORR_a_gene_pair' method.
    Return dataframe summarizing the results (swithced and not switched). If 'mask_zero_with_nan' is True, replace near-zero values with np.nan '''
    s_gene_pair = DEVIATION_NORMAL_CORR_a_gene_pair( gene_1 = gene_1, gene_2 = gene_2, df = df, n_std_for_outliers = n_std_for_outliers )
    s_gene_pair_switched = DEVIATION_NORMAL_CORR_a_gene_pair( gene_1 = gene_2, gene_2 = gene_1, df = df, n_std_for_outliers = n_std_for_outliers )
    df = pd.DataFrame( [ s_gene_pair, s_gene_pair_switched ], index = [ 'Initial', 'Switched' ] ).T
    df[ 'x' ] = df.Initial * df.Switched
    df[ '+' ] = df.Initial + df.Switched
    df[ '-' ] = df.Initial - df.Switched
    if mask_zero_with_nan : # if 'mask_zero_with_nan' is True, convert near zero values to np.nan
        df.values[ np.abs( df.values ) < 1e-8 ] = np.nan
    return df


# In[ ]:


def DEVIATION_NORMAL_CORR_a_gene_pair( gene_1, gene_2, df = None, n_std_for_outliers = 3.5 ) :
    ''' calculates 14 metrics for projection-based linear regression analysis for a given pair of genes on 'df' data. '''
    if df is None :
        df = df_proteome_unshared_dropna
    gene_id_1, gene_id_2 = Gene_2_Gene_ID( gene_1 ), Gene_2_Gene_ID( gene_2 )
    T_data_gene_1, T_data_gene_2 = GET_non_NaN_without_outliers( df[ T_samples ].loc[ gene_id_1 ], df[ T_samples ].loc[ gene_id_2 ], n_std_for_outliers = n_std_for_outliers )
    N_data_gene_1, N_data_gene_2 = GET_non_NaN_without_outliers( df[ N_samples ].loc[ gene_id_1 ], df[ N_samples ].loc[ gene_id_2 ], n_std_for_outliers = n_std_for_outliers )
    T_corr = NUMPY_spearman_correl( T_data_gene_1, T_data_gene_2 ) # calculate Spearman's rank correlation before projection
    N_corr = NUMPY_spearman_correl( N_data_gene_1, N_data_gene_2 )
    T_data_gene_1_2 = np.vstack( ( T_data_gene_1, T_data_gene_2 ) )
    N_data_gene_1_2 = np.vstack( ( N_data_gene_1, N_data_gene_2 ) )
    slope, intercept = LINREGRESS_Projection_Based_Least_Square_Linear_Regression( N_data_gene_1, N_data_gene_2, N_data_gene_1_2 )
    an_intercept_vector[ 1 ] = intercept # define vector of position of intercept
    R = np.array( [ [ 1, slope ], [ - slope, 1 ] ] ) / np.sqrt( 1 + slope ** 2 ) # calculation rotation matrix that rotate linear regression line to x-axis
    T_data_gene_1_2_projected = np.matmul( R, T_data_gene_1_2 - an_intercept_vector ) # subtract intercept_vector and perform projection to the linear regression line using the rotation matrix.
    N_data_gene_1_2_projected = np.matmul( R, N_data_gene_1_2 - an_intercept_vector )
    T_rej_mean, T_rej_std = NUMPY_GET_mean_std( T_data_gene_1_2_projected[ 1 ] ) # Than, calculating mean and standard deviation of rejections of the points 
    N_rej_mean, N_rej_std = NUMPY_GET_mean_std( N_data_gene_1_2_projected[ 1 ] )
    T_proj_mean, T_proj_std = NUMPY_GET_mean_std( T_data_gene_1_2_projected[ 0 ] ) # Than, calculating mean and standard deviation of rejections of the points 
    N_proj_mean, N_proj_std = NUMPY_GET_mean_std( N_data_gene_1_2_projected[ 0 ] )
    T_dev_corr = NUMPY_spearman_correl( T_data_gene_1_2_projected[ 0 ], T_data_gene_1_2_projected[ 1 ] ) # calculate Spearman's rank correlation between projections and rejections of points
    N_dev_corr = NUMPY_spearman_correl( N_data_gene_1_2_projected[ 0 ], N_data_gene_1_2_projected[ 1 ] )
    T_data_gene_1_2_projected[ 0 ].sort( ), N_data_gene_1_2_projected[ 0 ].sort( ) # sort the projection data in-place to retrive max and min of projection (this is one of the fastest way)
    T_proj_min, T_proj_max = T_data_gene_1_2_projected[ 0, 0 ], T_data_gene_1_2_projected[ 0, - 1 ]
    N_proj_min, N_proj_max = N_data_gene_1_2_projected[ 0, 0 ], N_data_gene_1_2_projected[ 0, - 1 ]
    origin_rej = intercept / np.sqrt( slope ** 2 + 1 )
    origin_proj = intercept * slope / np.sqrt( slope ** 2 + 1 )
    return pd.Series( dict( slope = slope, intercept = intercept, T_rej_mean = T_rej_mean, T_rej_std = T_rej_std, N_rej_mean = N_rej_mean, N_rej_std = N_rej_std, T_proj_mean = T_proj_mean, T_proj_std = T_proj_std, N_proj_mean = N_proj_mean, N_proj_std = N_proj_std, T_dev_corr = T_dev_corr, N_dev_corr = N_dev_corr, T_corr = T_corr, N_corr = N_corr, T_proj_min = T_proj_min, T_proj_max = T_proj_max, N_proj_min = N_proj_min, N_proj_max = N_proj_max, origin_rej = origin_rej, origin_proj = origin_proj ) )


# In[ ]:


def LINREGRESS_OPTIMIZE_search_slope_with_minumun_std_of_rejections( x, y, x_y_centered, num_probing = 5, num_repeat = 3 ) :
    ''' by using double for loop, search a slope that gives minimal standard deviation of rejections for a given data ('x_y_centered').
    in a search space, 'num_probing' number of search points were searched, and magnify the search space, and repeat the search for 'num_repeat' times
    about 3 times faster than Scipy.optimum.fmin function, though this function is less accurate than the function '''
    initial_slopes = NUMPY_least_square_linear_regress( x, y )[ 0 ], 1 / NUMPY_least_square_linear_regress( y, x )[ 0 ] # retrive initial slopes using least square linear regression methods
    lower_slope, upper_slope = sorted( initial_slopes ) # retrive lower and upper slope for setting points in search space using np.arange
    stds = np.zeros( num_probing + 2 ) # create an numpy array that will store search results (outputs of a function whose minimum will be searched)
    stds[ 0 ], stds[ -1 ] = LINREGRESS_OPTIMIZE_calculate_std_of_rejections( lower_slope, x_y_centered ), LINREGRESS_OPTIMIZE_calculate_std_of_rejections( upper_slope, x_y_centered ) # store std of slopes of lower and upper boundaries of the search space
    for index_repeat in range( num_repeat ) : # repeat search with more dense search points in smaller sub-search space
        step = ( upper_slope - lower_slope ) / ( num_probing + 1 ) # calculate step between search points
        slopes = np.arange( lower_slope, upper_slope + 1e-5, step ) # build a list of search points 
        for index_probing in range( 1, num_probing + 1 ) : # calculate output of the function (std of rejections at a given slope) for every search point (a slope in the search space)
            stds[ index_probing ] = LINREGRESS_OPTIMIZE_calculate_std_of_rejections( slopes[ index_probing ], x_y_centered )
        two_lowest_stds = stds.argsort( )[ : 2 ] # retrive indices of two lowest function outputs (std of rejections), assuming the two indices are of adjacent points and there is only one local minimum)
        slopes = slopes[ two_lowest_stds ] # retrive slopes of two lowest std of rejections
        stds[ 0 ], stds[ -1 ] = stds[ two_lowest_stds ][ slopes.argsort( ) ] # retrive std of rejections of lower_slope and upper_slope
        lower_slope, upper_slope = sorted( slopes ) # retrive lower and upper slope for setting search points for the next round of search
    return ( lower_slope + upper_slope ) / 2 # once all repeats of search ended, return the average of slopes with two lowest function outputs


# In[ ]:


def LINREGRESS_calculate_std_rejections_and_projections( data_1, data_2, slope, intercept ) :
    data_1_2 = np.vstack( ( data_1, data_2 ) )
    intercept_vector = np.array( [ [ 0 ], [ intercept ] ] ) # define vector of position of intercept
    R = np.array( [ [ 1, slope ], [ - slope, 1 ] ] ) / np.sqrt( 1 + slope ** 2 ) # calculation rotation matrix that rotate linear regression line to x-axis
    data_1_2_projected = np.matmul( R, data_1_2 - intercept_vector ) # subtract intercept_vector and perform projection to the linear regression line using the rotation matrix.
    return data_1_2_projected[ 1 ].std( ), data_1_2_projected[ 0 ].std( )


# In[ ]:


def Plot_vectors( data, label, alpha = 0.5 ) :
    plt.plot( data[ 0 ], data[ 1 ], 'o', alpha = alpha, label = label )
    plt.plot( np.arange( - 2, 2, 0.1 ), [ 0 ] * len( np.arange( - 2, 2, 0.1 ) )  , 'k' )
    plt.plot( [ 0 ] * len( np.arange( - 2, 2, 0.1 ) )  , np.arange( - 2, 2, 0.1 ), 'k' )    


# In[ ]:


def LINREGRESS_plot_slope_intercept( slope, intercept, data_1 = None, data_2 = None, x_start = 0, x_end = 2, alpha = 0.7, swithced_xy = False, label = '', calculate_std_rejections = False ) :
    x = np.arange( x_start, x_end, 0.03 )
    y = x * slope + intercept
    if calculate_std_rejections :
        rej_std, proj_std = LINREGRESS_calculate_std_rejections_and_projections( data_1, data_2, slope, intercept )
        label += '\nrejection std = {rej_std}'.format( rej_std = round( rej_std, 3 ) )
    else :
        label = 'a line of slope {slope}\nand intercept {intercept}'.format( slope = slope, intercept = intercept )
    if swithced_xy :
        plt.plot( y, x, '--', alpha = alpha, label = label )        
    else :
        plt.plot( x, y, '--', alpha = alpha, label = label )


# #### Solving a Cubic equation

# In[ ]:


"""def NUMPY_solve_cubic_equation( a, b, c, d ) :
    ''' solve cubic equation of a*x^3 + b*x^2 + c*x + d = 0. return a three solutions '''
    Q = ( 3 * a * c - b ** 2 ) / ( 9 * ( a ** 2 ) )
    R = ( 9 * a * b * c - 27 * ( a ** 2 ) * d - 2 * ( b ** 3 ) ) / ( 54 * ( a ** 3 ) ) 
    D = Q ** 3 + R ** 2
    if D < 1e-10 : # if D is Zero, since S and T is equal to np.cbrt( R ), calculate S-T and S+T directly to not use complex conjugate solution of np.power( negative real number + 0J, 1/3 ). np.power does not return real solution if input type is complex number type 
        S_plus_T = 2 * np.cbrt( R )
        S_minus_T = 0
    else : 
        sqrt_D = np.sqrt( D + 0J ) # convert D to complex number type to allow calculation of sqrt of negative D 
        S = NUMPY_CUBIC_FUNCTION_cuberoot_maximum_realpart_ratio( R + sqrt_D )
        T = NUMPY_CUBIC_FUNCTION_cuberoot_maximum_realpart_ratio( R - sqrt_D )
        S_minus_T = S - T
        S_plus_T = S + T
    b_by_3a = b / ( 3 * a )
    x_1 = S_plus_T - b_by_3a
    S_minus_T__times__i_sqrt3_by_2 = S_minus_T * 1J * np.sqrt( 3 ) / 2
    minus_b_by_3a__minus_S_plus_T__by_2 = - b_by_3a - S_plus_T / 2 
    x_2 = minus_b_by_3a__minus_S_plus_T__by_2 + S_minus_T__times__i_sqrt3_by_2
    x_3 = minus_b_by_3a__minus_S_plus_T__by_2 - S_minus_T__times__i_sqrt3_by_2
    return x_1, x_2, x_3, Q, R, S, T, D"""


# In[ ]:


"""def NUMPY_solve_cubic_equation_no_square_term( a, c, d ) :
    ''' solve cubic equation of a*x^3 + c*x + d = 0. return one real solution '''
    p_by_2 = c / a / 2
    q_by_3 = d / a / 3
    value_in_sqrt = p_by_2 ** 2 + q_by_3 ** 3
    if value_in_sqrt >= 0 :
        sqrt_value = np.sqrt( value_in_sqrt )
        x = np.cbrt( - p_by_2 + sqrt_value ) + np.cbrt( - p_by_2 - sqrt_value )
    else : 
        sqrt_value = np.sqrt( value_in_sqrt + 0J ) # calculate square root of negative number by converting the float datatype to complex datatype
        one_cuberoot = NUMPY_CUBIC_FUNCTION_cuberoot_maximum_realpart( - p_by_2 + sqrt_value ) 
        x = one_cuberoot + np.conj( one_cuberoot )
    return x, p_by_2, q_by_3, value_in_sqrt"""


# In[ ]:


"""def cuberoot( z ):
    x = z.real
    y = z.imag
    mag = abs(z)
    arg = math.atan2( y, x )
    resMag = mag**(1./3)
    resArg = [ ( arg + 2 * math.pi * n )/3. for n in range(1,4) ]
    return [  resMag*(math.cos(a) + math.sin(a)*1j) for a in resArg ]"""


# In[ ]:


"""def NUMPY_CUBIC_FUNCTION_cuberoot_maximum_realpart( z ) :
    ''' return a solution of cubic_root of complex number with maximal real part (for -0.2+1J and 0.1+1J, return the latter) '''
    x = z.real
    y = z.imag
    mag = abs(z)
    arg = math.atan2( y, x )
    resMag = mag**(1./3)
    resArg = [ ( arg + 2 * math.pi * n )/3. for n in range(1,4) ]
    real_part = np.array( [ math.cos( a ) for a in resArg ], dtype = float )
    return [  resMag*(math.cos(a) + math.sin(a)*1j) for a in resArg ][ real_part.argmax( ) ]"""


# In[ ]:


"""def NUMPY_CUBIC_FUNCTION_cuberoot_maximum_realpart_ratio( z ) :
    ''' return a solution of cubic_root of complex number with maximal real part contribution (for -0.2+1J and 0.1+1J, return the former) '''
    x = z.real
    y = z.imag
    mag = abs(z)
    arg = math.atan2( y, x )
    resMag = mag**(1./3)
    resArg = [ ( arg + 2 * math.pi * n )/3. for n in range(1,4) ]
    real_part_ratio = np.abs( np.array( [ math.cos( a ) for a in resArg ], dtype = float ) )
    return [  resMag*(math.cos(a) + math.sin(a)*1j) for a in resArg ][ real_part_ratio.argmax( ) ]"""


# ### Functions for Analysis of TCGA Data

# In[ ]:


def TCGA_CALCULATE_KaplanMeierCurve_log_rank_p_value( data_1, data_2, alpha = 0.05 ) :
    ''' For given two numpy arrays of survival data, calculate log rank p_value and return p_value, data_1_mean, data_2_mean, data_2_mean_by_data_1_mean  '''
    p_value = logrank_test( data_1, data_2, alpha = alpha ).p_value
    data_1_mean, data_2_mean = data_1.mean( ), data_2.mean( )
    return p_value, data_1_mean, data_2_mean, data_2_mean / data_1_mean


# In[ ]:


def TCGA_PLOT_KaplanMeierCurve( data_1, data_2, label_1 = 'data_1', label_2 = 'data_2', title = '', data_label = 'Days to Death (Days)', alpha = 0.05,  value_invalid = 5000, graph_folder = None, save_fig = False, show_fig = True ) :
    ''' For given two numpy arrays of survival data, data_1 and data_2, plot KaplanMeierCurve. Calculate p_value using 'TCGA_CALCULATE_KaplanMeierCurve_log_rank_p_value'
    alpha : type_1 error rate. value_invalid : used for setting x_limit  '''
    data_all = np.array( list( data_1 ) + list( data_2 ) ) # combine two data
    data_valid_entries = data_all[ data_all != value_invalid ] # retrive data of only valid entries (that is not 9999)
    x_lim = data_valid_entries.max( ) * 1.05 # set x_limit for the plot
    p_value, data_1_mean, data_2_mean, data_2_mean_by_data_1_mean = TCGA_CALCULATE_KaplanMeierCurve_log_rank_p_value( data_1, data_2, alpha = alpha ) # calculate p_value and ratio between data_2_mean and data_1_mean
    fig, ax = plt.subplots( 1, 1 )
    ax = KaplanMeierFitter( ).fit( data_1, label = label_1 ).plot( ax = ax )
    ax = KaplanMeierFitter( ).fit( data_2, label = label_2 ).plot( ax = ax )
    title = "{dataset_name} Kaplan-Meier Curves\n{title}\nLog Rank p-value= {p_value}".format( dataset_name = dataset_name, title = title, p_value = '%.2e' % p_value ) 
    ax.set_title( title )
    ax.set_xlabel( data_label )
    ax.set_xlim( left = 0, right = x_lim )
    ax.set_ylim( bottom = 0, top = 1.025 )
    if save_fig :
        plt.savefig( graph_folder + to_window_path_compatible_str( title ) + '.png', dpi = 200 )
        if not show_fig :
            plt.close( )
    else :
        return fig, ax


# In[ ]:


def TCGA_CALCULATE_KaplanMeierCurve_log_rank_p_value_of_all_genes_in_dataframe( df, Gene_Set = None, df_clinical = None, cut_ratio = 0.25 ) :
    ''' For each genes in df (or a subset of df if a 'Gene_Set' is given), calculate log rang p_values for data of 'days_to_death' column in df_clinical (default : df_patients) between
    two groups of patients, one group with upper n percent and another group with lower n percent of the data given by df '''
    if df_clinical is None : # set default df_clinical
        df_clinical = df_patients
    if Gene_Set is not None : # if 'Gene_Set' is given, subset df with a given Gene_Set
        df = PANDAS_Subset( df, Gene_Set )
    df = df[ T_samples ] # select only Tumor samples
    str_percent = str( int( round( cut_ratio * 100, 0 ) ) )
    data_1_name, data_2_name = 'lower_{}_percent'.format( str_percent ), 'upper_{}_percent'.format( str_percent ) # set data_1_name and data_2_name according to 'cut_ratio'
    Patient_IDs, _ = TCGA_SAMPLE_GET_List_Patients_IDs__List_Sample_Types__from__List_Sample_IDs( df.columns.values ) # retrive Patient_IDs from column names of a given df
    Gene_IDs, data = df.index.values, df.values # prepare list of Gene_IDs and data for iterations
    n_patients = len( Patient_IDs ) # retruve the number of patients
    dict_Gene_ID__analysis_result = dict( ) # an empty dictionary that will store analysis result
    for Gene_ID, data_of_a_gene in zip( Gene_IDs, data ) :
        data_of_a_gene__sorted = np.sort( data_of_a_gene ) # sort the data to retrive threshold values for lower n % and upper n % 
        thres_the_lower_portion = data_of_a_gene__sorted[ int( n_patients * cut_ratio ) ]
        thres_the_upper_portion = data_of_a_gene__sorted[ - int( n_patients * cut_ratio ) ]
        survival_data__low = df_clinical.loc[ Patient_IDs[ data_of_a_gene < thres_the_lower_portion ] ].days_to_death.values # retrive survival data of two groups of patients belonging to lower n_% and upper n_% of a given data of a gene
        survival_data__up = df_clinical.loc[ Patient_IDs[ data_of_a_gene > thres_the_upper_portion ] ].days_to_death.values
        p_value, data_1_mean, data_2_mean, data_2_mean_by_data_1_mean = TCGA_CALCULATE_KaplanMeierCurve_log_rank_p_value( survival_data__low, survival_data__up ) # calculate log_rank p_value
        dict_Gene_ID__analysis_result[ Gene_ID ] = dict( p_value = p_value, data_1_mean = data_1_mean, data_2_mean = data_2_mean, data_2_mean_by_data_1_mean = data_2_mean_by_data_1_mean )
    result = pd.DataFrame( dict_Gene_ID__analysis_result ).T.sort_values( 'p_value' ) # convert dictionary result into DataFrame and sort the result with p_values
    result = result.rename( columns = dict( data_1_mean = data_1_name + '__mean', data_2_mean = data_2_name + '__mean', data_2_mean_by_data_1_mean = 'Ratio of ' + data_2_name + ' to ' + data_1_name ) ) # rename the column names of the result dataframe
    return result # return the result DataFrame


# In[ ]:


def TCGA_PLOT_KaplanMeierCurve_for_a_gene( Gene, df = None, df_clinical = None, cut_ratio = 0.25, graph_folder = None, save_fig = False, show_fig = False ) :
    ''' For a given gene, draw a KaplanMeierCurve between two groups of patients in df_clinical (default : df_patients), one group with upper n percent and another group with lower n percent of the data given by df
    (default : df_rna). 'show_fig' : if both 'save_fig' and 'show_fig' are True, show the plot, and if 'save_fig' is False, always show the plot '''
    if df_clinical is None : # set default df_clinical
        df_clinical = df_patients
    if df is None :
        df = df_rna
    Gene_ID = Gene_2_Gene_ID( Gene ) # retrive Gene_ID from a given Gene
    if Gene_ID not in df.index.values :
        return -1
    Gene_Name_Symbol = dict_ID__Gene_Name_Symbol[ Gene_ID ]
    df = df[ T_samples ] # select only Tumor samples
    str_percent = str( int( round( cut_ratio * 100, 0 ) ) )
    data_1_name, data_2_name = 'low_{}_percent'.format( str_percent ), 'up_{}_percent'.format( str_percent ) # set data_1_name and data_2_name according to 'cut_ratio'
    Patient_IDs, _ = TCGA_SAMPLE_GET_List_Patients_IDs__List_Sample_Types__from__List_Sample_IDs( df.columns.values ) # retrive Patient_IDs from column names of a given df
    n_patients = len( Patient_IDs ) # retruve the number of patients
    data_of_a_gene = df.loc[ Gene_ID ].values
    data_of_a_gene__sorted = np.sort( data_of_a_gene ) # sort the data to retrive threshold values for lower n % and upper n % 
    thres_the_lower_portion = data_of_a_gene__sorted[ int( n_patients * cut_ratio ) ]
    thres_the_upper_portion = data_of_a_gene__sorted[ - int( n_patients * cut_ratio ) ]
    survival_data__low = df_clinical.loc[ Patient_IDs[ data_of_a_gene < thres_the_lower_portion ] ].days_to_death.values # retrive survival data of two groups of patients belonging to lower n_% and upper n_% of a given data of a gene
    survival_data__up = df_clinical.loc[ Patient_IDs[ data_of_a_gene > thres_the_upper_portion ] ].days_to_death.values
    return TCGA_PLOT_KaplanMeierCurve( survival_data__low, survival_data__up, label_1 = data_1_name, label_2 = data_2_name, title = Gene_Name_Symbol, graph_folder = graph_folder, save_fig = save_fig, show_fig = show_fig ) # plotted a Kaplan Meier Curve, and return fig and ax if save_fig is False


# ### CCLE Preprocessing functions

# In[ ]:


def Broad_CCLE_Cell_Line_Mapper( l_cl ) :
    ''' Map a given list of cell lines to CCLE Cell Line IDs. '''
    l_cl = list( cl.replace( '-', '' ).replace( ' ', '' ).upper( ) for cl in l_cl )
    for cl in l_cl : # map given list of cell lines to CCLE Cell line ids
        if cl not in df_ccle_meta_cl.index.values : # for cell line names that is not CCLE cell line ID
            l_cl_search = CCLE_Search_Cell_Line( cl + '_' )
            if len( l_cl_search ) > 1 :
                print( '[Cell Line Error] More than one cell line was searched :', l_cl_search, ' (Please refine your query)' )
            if len( l_cl_search ) == 0 :
                print( '[Cell Line Error] no cell line was found :', cl )
    return list( cl if cl in df_ccle_meta_cl.index.values else CCLE_Search_Cell_Line( cl + '_' )[ 0 ] for cl in l_cl )


# In[ ]:


def Broad_CCLE_Preprocessing_Meta__Cell_line_names( cl_name ) :
    cl_name_ccle = cl_name.replace( '-', '' ).replace( ':', '' ).split( '; ' )[ 0 ].upper( ).split( '(' )[ 0 ].strip( ).replace( ' ', '' )
    if cl_name_ccle in df_ccle_meta_cl.CL_Name.values :
        return cl_name_ccle 
    else :
        return np.nan


# In[ ]:


def Broad_CCLE_Preprocessing_map_unmapped_cell_lines( l_cl, print_result = True ) :
    ''' Return a dictionary that can rename cell_line_id to match the cell_line_id in df_ccle_meta_cl '''
    l_all_cl = df_ccle_meta_cl.index.values
    l_unmapped_cl = list( set( l_cl ).difference( l_all_cl ) )
    l_unmapped_cl_names = list( cl.split( '_' )[ 0 ] for cl in l_unmapped_cl )
    l_n_matched = list( np.array( list( True if cl_name == cl.split( '_' )[ 0 ] else False for cl in l_all_cl ) ).sum( ) for cl_name in l_unmapped_cl_names ) # number of ccle_cl that have exactly the same 
    df_matching_result = pd.DataFrame( dict( Unmapped_ID = l_unmapped_cl, Cell_Line = l_unmapped_cl_names, n_matched = l_n_matched ) ).set_index( 'Unmapped_ID' )
    if print_result :
        display( df_matching_result )
    df_matching_result = df_matching_result[ df_matching_result.n_matched == 1 ] # filter out invalid cell_line ids
    return dict( ( cl_id_unmaped, list( cl for cl in l_all_cl if cl_name == cl.split( '_' )[ 0 ] )[ 0 ] ) for cl_name, cl_id_unmaped in zip( df_matching_result.Cell_Line.values, df_matching_result.index.values ) )


# In[ ]:


def Broad_CCLE_Process_gene_dependency_data( df, fill_value_for_negative_values = 0.01 ) :
    ''' Process gene dependency data to make them positive before Log2FC calculation '''
    df = deepcopy( df )
    data = df.values
    data[ data > 0 ] = - fill_value_for_negative_values
    return pd.DataFrame( - data, columns = df.columns, index = df.index )


# ### CCLE Functions for Data Exploration

# In[ ]:


def CCLE_EXPORATORY_Add_snucellbank_availability( df, inplace = False ) :
    ''' Add SNU Cell Bank availability to a given dataframe where indices are CCLE cell line IDs  '''
    if not inplace :
        df = deepcopy( df )
    l_snu_cl = list( set( df_snucellbank_meta_cl.CCLE_ID.dropna( ) ) )
    df[ 'Available_in_SNU_Cell_Bank' ] = list( True if cl in l_snu_cl else False for cl in df.index.values )
    return df


# In[ ]:


def CCLE_EXPLORATORY_CL_Genes_Mutation_Table( l_genes = None, l_cl = None ) :
    ''' Return dataframe containing mutation status ( 'Protein_Change' column ) of a given list of cell lines. Cell line search result wil be displayed if there is error during searching
    Default cell liens are all the cell lines in ccle mutation database.
    TCGA Hotspot mutations are annotated with '(Hotspot)' mark at the end of mutation '''
    l_cl = df_mut.Tumor_Sample_Barcode.values if l_cl is None else Broad_CCLE_Cell_Line_Mapper( l_cl )
    l_genes = df_gdsc_cn.index.values if l_genes is None else l_genes # default genes are genes available in the copy number data
    l_geneids = List_Gene__2__List_Gene_ID( l_genes )
    dict_res = dict( )
    for geneid in l_geneids :
        df_search_res = CCLE_MUTATION_Retrive_a_gene_of_cell_lines( geneid, return_compact_df = False ).sort_values( 'isTCGAhotspot', ascending = False ).drop_duplicates( subset = [ 'Tumor_Sample_Barcode' ] ).set_index( 'Tumor_Sample_Barcode' ).astype( object ) # prioritize mutation data using TCGA_hotspot_mutation_status when more then protein_change data values are available
        arr_mut, arr_mask_hotspot, indices = df_search_res.Genome_Change.values, df_search_res.isTCGAhotspot.values.astype( bool ), df_search_res.index.values # retrive Genome Change values
        arr_mut[ ~ pd.isnull( df_search_res.cDNA_Change ).values ] = df_search_res.cDNA_Change.dropna( ).values # update cDNA Change values
        arr_mut[ ~ pd.isnull( df_search_res.Protein_Change ).values ] = df_search_res.Protein_Change.dropna( ).values # update Genome Change values
        arr_mut[ arr_mask_hotspot ] = arr_mut[ arr_mask_hotspot ] + ' (Hotspot)' # annotate hotspot mutation
        dict_res[ geneid ] = pd.Series( arr_mut, index = indices ) 
    df = PD_Subset( pd.DataFrame( dict_res ), l_cl ).fillna( 'WT' ).T
    for cl in set( l_cl ).difference( df.columns.values ) : # fill cell_lines with no mutations with 'WT'
        df[ cl ] = 'WT'
    df = df.T
    df_subset = PD_Subset( df, index = df_gdsc_cn.columns.values, columns = df_gdsc_cn.index.values ) # add df_gdsc_cn copy number data
    if df_subset.shape[ 1 ] != df.shape[ 1 ] :
        print( "[MAPPING ERROR] Following genes not present in 'df_gdsc_cn' dataset", List_Gene_ID__2__List_Gene_Symbol( set( df.columns.values ).difference( df_gdsc_cn.index.values ) ) )
    df_subset.values[ ( PD_Subset( df_gdsc_cn, index = l_geneids, columns = df_subset.index.values, preserve_order_in_df = False ) == 0 ).T & ( df_subset == 'WT' ) ] = 'DEL' # update copy number data (Deletion, Disruption, and Loss of Heterozygousy)
    df_subset.values[ PD_Subset( df_gdsc_cn_d, index = l_geneids, columns = df_subset.index.values, preserve_order_in_df = False ).T & ( df_subset == 'WT' ) ] = 'DIS'
    df_subset.values[ PD_Subset( df_gdsc_cn_loh, index = l_geneids, columns = df_subset.index.values, preserve_order_in_df = False ).T & ( df_subset == 'WT' ) ] = 'LOH'
    df = df_subset.append( df.loc[ list( set( df.index.values ).difference( df_subset.index.values ) ) ].replace( 'WT', 'WT?' ) ) # add back cell lines do not exist in df_gdsc_cn (they only have point mutation data) to the final processed result dataframe
    return PD_Add_gene_annotation( df.T, geneid_to_genesymbol = True ).T


# In[ ]:


def CCLE_Process_gene_dependency_data( df ) :
    df = deepcopy( df )
    data = df.values
    data[ data > 0 ] = - 0.01
    return pd.DataFrame( - data, columns = df.columns, index = df.index )


# In[ ]:


def CCLE_two_gene_plot( df, gene1, gene2 = None, df_2 = None, external_data = None, plot_kw = dict( ) ) :
    if external_data is not None :
        data_1 = PD_Locate_gene( df, gene1 )
        MPL_Scatter_Align_Two_Series( data_1, external_data, n_std_for_outliers = 0, annotate_xy_axis = ( True, True ), **plot_kw )   
    else :
        data_1 = PD_Locate_gene( df, gene1 )
        data_2 = PD_Locate_gene( df, gene2 ) if df_2 is None else PD_Locate_gene( df_2, gene2 )
        MPL_Scatter_Align_Two_Series( data_1, data_2, n_std_for_outliers = 0, annotate_xy_axis = ( True, True ), **plot_kw )


# In[ ]:


def CCLE_MUTATION_Retrive_a_cell_line( cell_line_name, return_compact_df = True ) :
    df = PD_Select( df_mut, Tumor_Sample_Barcode = cell_line_name )
    if return_compact_df :
        df = df[ [ 'Hugo_Symbol', 'Entrez_Gene_Id', 'Tumor_Sample_Barcode', 'Variant_Classification', 'Codon_Change', 'Protein_Change', 'isTCGAhotspot' ] ]
    return df


# In[ ]:


def CCLE_MUTATION_Retrive_a_gene_of_cell_lines( gene, l_cl = None, return_compact_df = True ) :
    df = PANDAS_Select( df_mut, Entrez_Gene_Id = Gene_2_Gene_ID( gene ) ) if l_cl is None else PANDAS_Select( df_mut, Entrez_Gene_Id = Gene_2_Gene_ID( gene ), Tumor_Sample_Barcode = l_cl ) # if valid list of cell_lines were given, subset df with given list of cell lines
    if return_compact_df :
        df = df[ [ 'Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification', 'Codon_Change', 'Protein_Change', 'isTCGAhotspot' ] ]
    return df


# In[ ]:


def CCLE_Search_Cell_Line( query ) :
    return Search_list_of_strings( df_ccle_meta_cl.index.values, query.upper( ) )


# In[ ]:


def CCLE_METABOLOMICS_Search_Metabolite( query ) :
    query = query.lower( )
    return sorted( list( metabol for metabol in df_metabol.index.values if query in metabol.lower( ) ) )


# In[ ]:


def CCLE_Search_Cell_Lines_df( df, query ) :
    query = query.upper( )
    return sorted( list( col for col in df.columns.values if query in col.upper( ) ) )


# In[ ]:


def CCLE_Search_Datasets_A_Cell_Line( cell_line_name ) :
    return list( data_name for df, data_name in zip( [ df_rna, df_rnai_gd, df_crispr_gd, df_metabol ], [ 'df_rna', 'df_rnai_gd', 'df_crispr_gd', 'df_metabol' ] ) if cell_line_name in df.columns.values )


# In[ ]:


def CCLE_Gene_Dependency_Score_Plot_Gene( Gene, data_type = 'crispr', save_fig = False, dict_cell_line_set = None, gene_for_rna = None, size_factor = 1, thres_n_cell_line = 5, show_legend = True ) :
    '''  Visualize Gene_Dependency data along with RNA data. bubble size indicates RNA amount (linear), and can be adjusted by 'size_factor'. Sort cell_line_set according to its median gene_dependency 
    value, while excluding cell_line_set if a subset contain less than 'thres_n_cell_line' entries, 'gene_for_rna' : Gene symbol for gene used for retriving RNA-Seq data. by default, the given gene for retriving GD Data   '''
    if dict_cell_line_set is None :
        dict_cell_line_set = dict_ccle__cell_line_sets
    df = df_crispr_gd if data_type.lower( ) == 'crispr' else df_rnai_gd # set df according to 'data_type'
    Gene_ID, Gene_Symbol = GET_Gene_ID_and_Symbol_from_input_gene( Gene )
    geneid_for_rna = Gene_ID if gene_for_rna is None else Gene_2_Gene_ID( gene_for_rna )
    if Gene_ID == -1 or Gene_ID not in df.index.values : # if gene is invalid or not exist in the data, return an error value
        return -1
    df = pd.DataFrame( dict( Gene_Dependency = df.loc[ Gene_ID ].dropna( ), RNA = df_rna.loc[ geneid_for_rna ] ) ).dropna( )
    size_factor = 1 / np.median( df.RNA.values ) * 20 * size_factor
    x_range = np.arange( - len( df ) * 0.1, len( df ) * 1.03 ) # plot x axis
    fig, ax = plt.subplots( 1, figsize = ( 12, 5 ) )
    ax.plot( x_range, np.zeros_like( x_range ), color = 'k', lw = 1 )
    ax.plot( x_range, np.full_like( x_range, -1 ), color = 'r', lw = 1 )
    s_cell_line_set_median = pd.Series( dict( ( cell_line_set_name, np.median( PANDAS_Subset( df, cell_line_set ).Gene_Dependency ) ) for cell_line_set_name, cell_line_set in dict_cell_line_set.items( ) if len( PANDAS_Subset( df, cell_line_set ) ) > thres_n_cell_line ) ) # retrive median viability score for each cell_line_set. exclude cell_line_set if a subset contain less than 'thres_n_cell_line' entries 
    x_position = 0
    for cell_line_set_name in s_cell_line_set_median.dropna( ).sort_values( ).index.values :
        cell_line_set = dict_cell_line_set[ cell_line_set_name ]
        df_subset = PANDAS_Subset( df, cell_line_set ).sort_values( 'Gene_Dependency' )
        ax.scatter( x = np.arange( len( df_subset ) ) + x_position, y = df_subset.Gene_Dependency, s = df_subset.RNA.values * size_factor, label = cell_line_set_name, alpha = 0.5 )
        x_position += len( df_subset )
    ax.set_ylim( - 1.5, 0.5 )
    ax.set_xlim( - len( df ) * 0.05, len( df ) * ( 1 + show_legend * 0.7 ) )
    ax.set_xticks( np.arange( 0, len( df ) + 100, 100 ) )
    data_description, data_method = ( 'CRISPR KO Viability (CERES Data)', 'Knock Out' ) if data_type == 'crispr' else ( 'RNAi KD Viability (DEMETER2 Data)', 'Knock Down' ) # set data_description and data_abbreviation according to 'data_type'
    MATPLOTLIB_basic_configuration( x_label = 'Ranked Cell Lines', y_label = 'Viability Score After {}'.format( data_method ), title = "'{}' {}\nCell Lines By Primary Site".format( Gene_Symbol, data_description ), save_fig = save_fig, show_legend = True )


# In[ ]:


def CCLE_CTRP_GD_Plot( cpd, gene, save_fig = False ) :
    gene_id, gene_symbol = GET_Gene_ID_and_Symbol_from_input_gene( gene )
    if 'df_ctrp_auc_aligned' not in globals( ) : # if aligned data is not available, align the two data
        globals( )[ 'df_ctrp_auc_aligned' ], globals( )[ 'df_crispr_gd_aligned' ] = PANDAS_Align_two( df_ctrp_auc, df_crispr_gd )
    x_range = np.arange( 0, 15.5, 0.1 )
    y_range = np.zeros_like( x_range )
    plt.plot( x_range, y_range, 'black', lw = 1 )
    plt.plot( x_range, y_range - 1, 'red', lw = 1 )
    plt.plot( df_ctrp_auc_aligned.loc[ cpd ].values, df_crispr_gd_aligned.loc[ gene_id ].values, 'o', alpha = 0.5 )
    MATPLOTLIB_basic_configuration( y_lim = ( -1.5, 0.5 ), x_label = "Area Under Curve (AUC) of '{}'".format( cpd ), y_label = "Gene Dependency Data of '{}'".format( gene_symbol ), title = "{} and {}\n(CTRP Drug Response Data, CRISPR KO Gene Dependency Data)".format( cpd, gene_symbol ), save_fig = save_fig, save_file_name = "{} and {}_(CTRP AUC_CRISPR GD)".format( cpd, gene_symbol ) )


# In[ ]:


"""def CCLE_CRISPR_Score_Plot_Gene( Gene_Symbol, save_fig = False ) :
    data = df_crispr_gd.loc[ dict_Symbol_ID_simple[ Gene_Symbol ] ].sort_values( ).dropna( ).values
    x_range = np.arange( len( data ) )
    plt.plot( x_range, np.zeros_like( x_range ), color = 'k', lw = 1 )
    plt.plot( data )
    plt.ylim( -1, 0.5 )
    MATPLOTLIB_basic_configuration( x_label = 'Ranked Cell Lines', y_label = 'Viability Score After KO', title = '{} CRISPR KO Viability in 558 Cell Lines'.format( Gene_Symbol ), save_fig = save_fig )"""


# In[ ]:


def CCLE_Scatter_RNA_Seq_and_Metabolite( Gene, Metabolite, query = None, log_scale = True ) :
    ''' 'query' = '?' : return list of tissues and histology subtypes '''
    df = pd.DataFrame( dict( RNA_Seq = df_rna.loc[ dict_Symbol_ID_simple[ Gene ] ], Metabolomics = df_metabol.loc[ Metabolite ] ) ).dropna( )
    if query == '?' : # if help is required, print out tissue types and histology subtypes
        print( 'Tissue Types : ' )
        display( sorted( set( list( index.split( '_' )[ 1 ].capitalize( ) for index in df_meta_ccle.index.values ) ) ) )
        print( 'Histology Subtypes : ' ) 
        display( df_meta_ccle.Hist_Subtype1.dropna( ).drop_duplicates( ).sort_values( ).values )
        return -1 
    if type( query ) is str : # if only tissue is given
        df = df.loc[ CCLE_Search_Cell_Lines_df( df.T, query ) ]
    elif query is not None : # if tissue and subtype query is given
        df_meta_ccle_matched_tissue = PANDAS_Subset_df( df_meta_ccle.dropna( subset = [ 'Hist_Subtype1' ] ), CCLE_Search_Cell_Lines_df( df.T, query[ 0 ] ) )
        df_subtype_searched = df_meta_ccle_matched_tissue[ list( True if query[ 1 ].lower( ) in subtype.lower( ) else False for subtype in df_meta_ccle_matched_tissue.Hist_Subtype1.values ) ]
        display( list( df_subtype_searched.Hist_Subtype1.drop_duplicates( ).sort_values( ).values )  ) # display searched histology subtypes
        df = PANDAS_Subset_df( df, df_subtype_searched.index.values )
    correl_p, p_p = stats.pearsonr( df.RNA_Seq.values, df.Metabolomics.values ) # calculate Pearson's correlation
    correl_s, p_s = stats.spearmanr( df.RNA_Seq.values, df.Metabolomics.values ) # calculate Spearman's correlation
    plt.plot( df.RNA_Seq.values, df.Metabolomics.values, '.' )
    if log_scale :
        scale = 'log'
    else :
        scale = 'linear'
    MATPLOTLIB_basic_configuration( x_scale = scale, y_scale = scale, x_label = 'RNA-Seq ({})'.format( Gene ), y_label = 'Metabolite Amount ({})'.format( Metabolite ), show_legend = False, title = 'RNA-Seq and Metabolomic Data', save_fig = False  )
    print( '\n(Correl. Pearson = {c_p}, Spearman = {c_s} \n P-value. Pearson = {p_p}, Spearman = {p_s})'.format(c_p = round( correl_p, 3 ), c_s = round( correl_s, 3 ), p_p = round( p_p, 3 ), p_s = round( p_s, 3 ) ) )


# ##### CCLE functions PCA and t-SNE

# In[ ]:


def CCLE_Transform_to_Standard_Score( df, l_cl = None, transform_log2 = False, transform_z_score = True, transform_gd_data = True ) :
    ''' Transform datavalues to standard scores for the use in PCA or t-SNE  '''
    if l_cl is not None :
        df = PD_Subset( df, l_cl, axis = 1 )
    if transform_gd_data :
        CCLE_Process_gene_dependency_data( df )
    if transform_log2 :
        df = np.log2( df.replace( 0, 1e-3 ) )
    if transform_z_score :
        df = STANDARD_SCORE_Z_score( df )
    return df


# ### IMPC_CRBN_KO functions

# In[ ]:


def IMPC_CRBN_KO__BIOTYPE_add_biotype_to_transcript_dataframe( df ) :
    df = deepcopy( df )
    cols = list( df.columns.values )
    if type( df.index.values[ 0 ] ) is tuple : # if the given dataframe is multiIndexed
        arr_tx_id = PANDAS_MULTIINDEX_get_indices_from_multiIndex_of_df( df, 1 )
    else :
        arr_tx_id = df.index.values
    arr_gene_biotype, arr_transcript_biotype = df_transcript_biotype.loc[ arr_tx_id ].values.T
    df[ 'gene_biotype' ] = arr_gene_biotype
    df[ 'transcript_biotype' ] = arr_transcript_biotype
    return df[ [ 'gene_biotype', 'transcript_biotype' ] + cols ]


# In[ ]:


def IMPC_CRBN_KO__COMBINE_Log2FC_results( df_result_liver, df_result_lung, df_result_heart, col_essential = [ 'Log2_Fold_Change', 'Difference', 'p_value', 'WT_average', 'CRBN_KO_average' ] ) :
    col_essential_all, col_essential_renamed_liver, col_essential_renamed_lung, col_essential_renamed_heart = list( col + '_all' for col in col_essential ), list( col + '_liver' for col in col_essential ), list( col + '_lung' for col in col_essential ), list( col + '_heart' for col in col_essential )
    df_result_liver, df_result_lung, df_result_heart = df_result_liver[ col_essential ], df_result_lung[ col_essential ], df_result_heart[ col_essential ]
    df_result_liver, df_result_lung, df_result_heart = df_result_liver.rename( columns = dict( ( col, col + '_liver' ) for col in col_essential ) ), df_result_lung.rename( columns = dict( ( col, col + '_lung' ) for col in col_essential ) ), df_result_heart.rename( columns = dict( ( col, col + '_heart' ) for col in col_essential ) )
    df_result_all = df_result_liver.join( df_result_lung, how = 'outer' ).join( df_result_heart, how = 'outer' )
    arr_cols = np.array( [ col_essential_renamed_liver, col_essential_renamed_lung, col_essential_renamed_heart ] ).T
    cols_summary = list( )
    for col_all, cols_tissues in zip( col_essential_all, arr_cols ) :
        if 'p_value' in col_all :
            df_result_all[ col_all + '_log_avg' ] = np.power( df_result_all[ cols_tissues ].prod( axis = 1 ), 1 / 3 ) # calculate geometric mean of p_values
            df_result_all[ col_all ] = df_result_all[ cols_tissues ].mean( axis = 1 ) # calculate linear mean of p_values
            cols_summary.append( col_all + '_log_avg' )
        else :
            df_result_all[ col_all ] = df_result_all[ cols_tissues ].mean( axis = 1 ) # calculate linear average for all columns other than p_values
        cols_summary.append( col_all )
    return df_result_all[ cols_summary + list( arr_cols.ravel( ) ) ]


# ### Functions for Nature_Ubi_Profiling 

# In[ ]:


def UBIQUITINATION_DATASET_barplot_a_protein( Gene, save_fig = False, df = None, data_use_log2 = False, CRBN_only = True ) :
    ''' Using df_ubi, draw barplots of data_values for a given gene '''
    if df is None :
        if CRBN_only :
            df = df_ubi[ ( cols_E1_E2 + cols_E1_E2_CRBN + cols_E1_E2_CRBN_Len + cols_E1_E2_CRBN_CSN )[ : : -1 ] ]
        else :
            df = df_ubi[ ( cols_E1_E2 + cols_E1_E2_CRBN + cols_E1_E2_CRBN_Len + cols_E1_E2_CRBN_CSN + cols_E1_E2_Cdt2 + cols_E1_E2_DDB2 )[ : : -1 ] ]
    Gene_ID_Symbol = GET_Gene_ID_and_Symbol_from_input_gene( Gene )
    if Gene_ID_Symbol == -1 :
        return -1
    Gene_ID, Gene_Symbol = Gene_ID_Symbol
    data = df.loc[ Gene_ID ]
    labels = data.index.values
    data_values = data.values
    if data_use_log2 :
        data_values = np.log2( data_values )
    plt.figure( figsize = ( 7, 5 ) )
    plt.barh( labels, data_values )
    MATPLOTLIB_basic_configuration( x_label = 'Ubiquitination Signal' + ' (Log2)' * data_use_log2, title = "Ubi. Signal of ( {} )".format( Gene_Symbol ) + ' by CRBN' * CRBN_only + "\n(2014 Nature Fischer Dataset)", save_fig = save_fig )


# In[ ]:


def UBIQUITINATION_DATASET_barplot_proteins_for_pub( Genes, background_method = 'background_divide', save_fig = False, adjust_aspect = 1, log_scale = False ) :
    ''' divide or subtract df_ubi data with its background data (E1_E2) according to 'background_method', and draw barplots using the data for a given list of genes for publication purposes '''
    Genes = LIST_intersection_with_set( List_Gene__2__List_Gene_ID( Genes ), df_ubi.index.values )
    print( len( Genes ) )
    if background_method == 'background_divide' : # divide df_ubi data with its background data (E1_E2)
        df = ( df_ubi.T / df_ubi[ cols_E1_E2 ].values.mean( axis = 1 ) ).T - 1
        aspect = 0.3 * adjust_aspect
    else : # subtract df_ubi data with its background data (E1_E2)
        df = ( df_ubi.T - df_ubi[ cols_E1_E2 ].values.mean( axis = 1 ) ).T
        aspect = 50 * adjust_aspect
    fig, arr_ax = plt.subplots( len( Genes ), sharex = True, figsize = ( 5, len( Genes ) ) )
    for Gene, ax in zip( Genes, arr_ax ) :
        Gene_ID_Symbol = GET_Gene_ID_and_Symbol_from_input_gene( Gene )
        if Gene_ID_Symbol == -1 :
            return -1
        Gene_ID, Gene_Symbol = Gene_ID_Symbol
        data = df.loc[ Gene_ID ]
        mean_E1_E2_CRBN, err_E1_E2_CRBN = NUMPY_GET_mean_std( data.loc[ cols_E1_E2_CRBN ].values )
        mean_E1_E2_CRBN_Len, err_E1_E2_CRBN_Len = NUMPY_GET_mean_std( data.loc[ cols_E1_E2_CRBN_Len ].values )
        labels = [ 'Lenalidomide + CUL4-CRBN', ' CUL4-CRBN']
        data_values = [ mean_E1_E2_CRBN_Len, mean_E1_E2_CRBN ]
        error_values = [ err_E1_E2_CRBN_Len, err_E1_E2_CRBN ]
        ax.barh( labels, data_values, xerr = error_values, capsize = 3 ) # draw bar plots with error values
        ax.set_title( Gene_Symbol )
        if log_scale :
            ax.set_xscale( 'log' )
        else :
            ax.set_aspect( aspect )
        ax.grid( False )
    plt.xlabel( 'Ubiquitination Signal Relative to Background' )
    title = "Ubiquitination Signals"
    fig.suptitle( title )
    if save_fig : # if save_fig is True, save a plot with title as its filename
        MATPLOTLIB_savefig( title = title + '_' + TIME_GET_timestamp( ) )


# In[ ]:


def CRBNPROT_Decompress_metadata( df = None, data_labels = None ) : 
    ''' Decompress metadata and Create metadata dataframe from columns of the given dataframe containing CRBN_Proteomic profiling dataset. Also see CRBNPROT_Compress_metadata
    Default data_labels : [ 'Value_Type', 'Data_Type', 'Experiment_System', 'Condition_A', 'Condition_B', 'Replicate', 'Data_Source' ] '''
    data_labels = [ 'Value_Type', 'Data_Type', 'Experiment_System', 'Condition_A', 'Condition_B', 'Replicate', 'Data_Source' ] if data_labels is None else data_labels # set default data_labels
    df = df_crbnprot if df is None else df # default df is df_crbnprot
    return pd.DataFrame( list( col.split( '__' ) for col in df_crbnprot.columns.values ), index = df_crbnprot.columns.values, columns = data_labels )


# In[ ]:


def CRBNPROT_Compress_metadata( df = None, data_labels = None ) : 
    ''' Compress metadata in the given metadata dataframe (default, df_crbnprot_meta_assays) according to data_labels and return list of labels containing compressed metadata.
    Default data_labels : [ 'Value_Type', 'Data_Type', 'Experiment_System', 'Condition_A', 'Condition_B', 'Replicate', 'Data_Source' ]. Also see CRBNPROT_Decompress_metadata '''
    data_labels = [ 'Value_Type', 'Data_Type', 'Experiment_System', 'Condition_A', 'Condition_B', 'Replicate', 'Data_Source' ] if data_labels is None else data_labels # set default data_labels
    df = df_crbnprot_meta_assays if df is None else df # default df is df_crbnprot_meta_assays
    df = PD_Subset( df, data_labels, axis = 1 ) # subset columns of df with data_labels 
    return list( '__'.join( values ) for values in df.values ) # return compressed metadata


# In[ ]:


def CRBNPROT_Locate_gene( gene, add_decompressed_metadata = True ) :
    ''' Retrive data of a given gene in the CRBN Proteome Profiling datasets '''
    s_data = PD_Locate_gene( df_crbnprot, gene )    
    if type( s_data ) is not pd.Series : 
        return -1
    geneid, genesymbol = GET_Gene_ID_and_Symbol_from_input_gene( gene )
    labels = s_data.index.values
    l_labels_log2fc, l_labels_p_values = Search_list_of_strings( labels, 'log2fc' ), Search_list_of_strings( labels, 'value' )
    df_meta = PD_Subset( PD_Subset( df_crbnprot_meta_assays, l_labels_log2fc ), [ 'Data_Type', 'Experiment_System', 'Condition_A', 'Condition_B' ], axis = 1 ) # subset metadata dataframe
    if add_decompressed_metadata :
        df = pd.DataFrame( dict( Log2_Fold_Change = s_data.loc[ l_labels_log2fc ].values, p_value = s_data.loc[ l_labels_p_values ].values ), index = l_labels_log2fc ).join( df_meta )
        df[ 'index' ] = CRBNPROT_Compress_metadata( df_meta )
        df = df.set_index( 'index' ) # set_index of the result df
    else :
        df = pd.DataFrame( dict( Log2_Fold_Change = s_data.loc[ l_labels_log2fc ].values, p_value = s_data.loc[ l_labels_p_values ].values ), index = CRBNPROT_Compress_metadata( df_meta ) )
    df.index.name = "Data for {}".format( genesymbol ) # set the name of the dataframe as the name of its index
    return df


# In[ ]:


def CRBNPROT_Retrive_list_of_log2fc_and_p_value_cols( df = None ) :
    ''' Return log2fc, p_value columns of a given df (default: df_crbnprot) '''
    cols = df_crbnprot.columns.values if df is None else df.columns.values
    return Search_list_of_strings( cols, 'log2fc' ), Search_list_of_strings( cols, 'value' ) # and 'Ubi_P' not in col and 'Author' not in col )


# ### Functions for online tools

# In[ ]:


def ONLINE_TOLLS_Enrichr( df, col = 'Log2_Fold_Change' ) :
    ''' To utilize public resource Enrichr, prepare text file (by ipynb output) for inpur of Enrichr  '''
    df = PD_Add_gene_annotation( df )
    for symbol, weight in zip( df.Approved_Symbol.values, BOKEH_transform_arr_values_into_a_range_from_0_to_1( df[ col ].values ) ) :
        print( "{},{}".format( symbol, weight ) )


# ### External Functions (from stack overflow)

# In[ ]:


def heatmap(data, row_labels, col_labels, show_colorbar = True, ax=None, cbar_kw={}, cbarlabel = None, color_limit = None, grid_color = 'w', grid_width = 3, fontsize_x = 10, fontsize_y = 10, **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    color_limit = ( low, high ) tuple
    """
    if not ax:
        ax = plt.gca()
    im = ax.imshow(data, **kwargs) # Plot the heatmap
    if show_colorbar : # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    if color_limit is None :
        pass
    else :
        im.set_clim( color_limit[ 0 ], color_limit[ 1 ] )
    if cbarlabel is not None and show_colorbar : # label color bar if name has been given
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom" )
    ax.set_xticks(np.arange(data.shape[1]))     # We want to show all ticks...
    ax.set_yticks(np.arange(data.shape[0]))
    ax.set_xticklabels( col_labels, fontsize = fontsize_x ) # ... and label them with the respective list entries.
    ax.set_yticklabels( row_labels, fontsize = fontsize_y )
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False) # Let the horizontal axes labeling appear on top.
    plt.setp( ax.get_xticklabels( ), rotation=-30, ha="right", rotation_mode="anchor") # Rotate the tick labels and set their alignment.
    for edge, spine in ax.spines.items(): # Turn spines off and create white grid.
        spine.set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color= grid_color, linestyle='-', linewidth = grid_width )
    ax.tick_params(which="minor", bottom=False, left=False)
    if show_colorbar :
        return im, cbar
    else :
        return im, None


# In[ ]:


def annotate_heatmap(im, data=None, valfmt="{x:.2f}", textcolors=["black", "white"], threshold=None, threshold_bothends=None, **textkw):
    """    A function to annotate a heatmap.
    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.
        threshold_bothends : a proportion of a diverging colormap according to which the colors from
                             textcolors are applied.
    Further arguments are passed on to the created text labels.    """
#     if not isinstance( data, (list, np.ndarray ) ) : # if data has not given, retrive data from the given image
    data_heatmap = im.get_array( ) # retrive data used to draw heatmap from the given image
    if threshold_bothends is None : # Normalize the threshold to the images color range.
        if threshold is not None:
            threshold = im.norm( threshold )
        else:
            threshold = im.norm( data_heatmap.max( ) )/2.
    else :
        threshold_low = im.norm( data_heatmap.max( ) ) * ( 1. - threshold_bothends )
        threshold_high = im.norm( data_heatmap.max( ) ) * threshold_bothends
    # Set default alignment to center, but allow it to be
    kw = dict( horizontalalignment = "center", verticalalignment = "center" ) # overwritten by textkw.
    kw.update( textkw )
    # Get the formatter in case a string is supplied
    if isinstance( valfmt, str ):
        valfmt = mpl.ticker.StrMethodFormatter( valfmt )
    # Loop over the data and create a `Text` for each "pixel".
    texts = [ ] # Change the text's color depending on the data.
    if threshold_bothends is None :
        for i in range( data.shape[ 0 ] ):
            for j in range( data.shape[ 1 ] ):
                kw.update(color = textcolors[ im.norm( data_heatmap[ i, j ] ) < threshold ])
                text = im.axes.text(j, i, valfmt( data[i, j], None ), **kw ) # annotate data
                texts.append( text )
    else : 
        for i in range( data.shape[ 0 ] ):
            for j in range( data.shape[ 1 ] ):
                kw.update( color = textcolors[ im.norm( data_heatmap[ i, j ] ) > threshold_high | im.norm( data_heatmap[ i, j ] ) < threshold_low ] )
                text = im.axes.text( j, i, valfmt(data[ i, j ], None ), **kw ) # annotate data
                texts.append( text )
    return texts


# In[ ]:


def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        # print ( chrom )
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


# ### Functions parsing files into DataFrames

# In[ ]:


def GTF_Parse_Attribute( attr ) :
    """
    # 2021-02-06 18:51:47 
    parse attribute string of a gtf file
    """
    dict_data = dict( )
    for e in attr.split( ';' ) :
        e = e.strip( )
        if len( e ) == 0 : 
            continue
        str_key, str_value = e.split( ' "' )
        if str_value[ -1 ] == '"' :
            str_value = str_value[ : -1 ]
        dict_data[ str_key ] = str_value
    return dict_data

def GTF_Read( dir_gtf, flag_gtf_gzipped = False, parse_attr = False ) :
    ''' 
    # 2021-02-06 18:50:59 
    Load gzipped or plain text GTF files into pandas DataFrame. the file's gzipped-status can be explicitly given by 'flag_gtf_gzipped' argument. 
    'parse_attr' : parse gtf attribute if set to True
    '''
    
    try :
        df = pd.read_csv( dir_gtf, sep = '\t', header  = None, low_memory = False ) # if there is no comments
    except pd.errors.ParserError : # if there is several comments, remove the comments and read the file
        if flag_gtf_gzipped is False and dir_gtf[ - 3 : ].lower( ) == '.gz' : flag_gtf_gzipped = True # automatically detect the file's gzipped status if it has not been explicitly given
        if flag_gtf_gzipped : # open VCF file based on its gzipped-status
            with gzip.open( dir_gtf, 'rb' ) as file :
                while True :
                    line = file.readline( ).decode( 'utf-8' )
                    if line[ 0 ] != '#' : break
                df = pd.read_csv( StringIO( line + '\n' + file.read( ).decode( 'utf-8' ) ), sep = '\t', header  = None, low_memory = False )
        else :
            with open( dir_gtf, 'r' ) as file :
                while True :
                    line = file.readline( )
                    if line[ 0 ] != '#' : break
                display( line )
                df = pd.read_csv( StringIO( line + file.read( ) ), sep = '\t', header  = None, low_memory = False )
    df.columns = [ 'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute' ]
    df = df.sort_values( [ 'seqname', 'start' ] ).reset_index( drop = True )
    if parse_attr :
        return df.join( pd.DataFrame( list( GTF_Parse_Attribute( attr ) for attr in df.attribute.values ) ) )
    return df


# In[ ]:


def PARSER_EMBL_Format( dir_file ) :
    ''' Parse EMBL format into DataFrame '''
    with open( dir_file, 'r' ) as file :
        l_record = file.read( ).split( '//\n' )

    l_dict_field_to_content = list( )
    for record in l_record[ : -1 ] :
        l_record_lines = record.split( '\n' )
        dict_field_to_content = dict( )
        for line in record.split( '\n' ) :
            str_field = 'Sequence' if line[ : 2 ] == '  ' else line[ : 2 ] # if field is empty, assign its field to 'Sequence' field
            str_content = line.split( '{}   '.format( line[ : 2 ] ) )[ 1 ] if '   ' in line else ''
            dict_field_to_content[ str_field ] = dict_field_to_content[ str_field ] + ' ' + str_content if str_field in dict_field_to_content else str_content
        l_dict_field_to_content.append( dict_field_to_content )
    df = pd.DataFrame( l_dict_field_to_content ).drop( columns = [ 'XX' ] )
    return df


# ### Functions for genomic data

# In[ ]:


def GTF_Calculate_overlapping_percentage( start1, end1, start2, end2 ) :
    ''' calculate and return overlapping percentages of two pairs of start and end coordinates  '''
    if start1 <= start2 :
        if end1 >= end2 :
            return ( end2 - start2 ) / ( end1 - start1 ) * 100, 100
        else :
            if end1 < start2 :
                return 0, 0
            else :
                overlap_length = end1 - start2
                return overlap_length / ( end1 - start1 ) * 100, overlap_length / ( end2 - start2 ) * 100
    else :
        if end2 >= end1 :
            return 100, ( end1 - start1 ) / ( end2 - start2 ) * 100
        else :
            if end2 < start1 :
                return 0, 0
            else :
                overlap_length = end2 - start1
                return overlap_length / ( end1 - start1 ) * 100, overlap_length / ( end2 - start2 ) * 100


# #### Functions for searching sequences

# In[ ]:


# from collections import defaultdict

def Substring_allow_mismatch( pattern, text ) :
    """ Return True if pattern exist in text allowing upto one mismatch.
    source : https://stackoverflow.com/questions/2420412/search-for-string-allowing-for-one-mismatch-in-any-location-of-the-string """
    m = len( pattern )
    S_table = defaultdict( int )
    for i, c in enumerate( pattern ) :
        S_table[c] |= 1 << i
    R0 = 0
    R1 = 0
    mask = 1 << ( m - 1 )
    for j, c in enumerate( text ) :
        S = S_table[c]
        shR0 = ( R0 << 1 ) | 1
        R0 = shR0 & S
        R1 = ( ( R1 << 1 ) | 1 ) & S | shR0
        if R0 & mask or R1 & mask : 
            return True
    return False


# ### Functions for interacting with OS

# In[ ]:


def OS_Currently_running_processes( ) :
    ''' Parse the output of a command line 'ps -ef' and parse it into a DataFrame, and return the DataFrame '''
    l_lines = os.popen( 'ps -ef' ).read( ).split( '\n' )[ 1 : -1 ]
    l_csv_lines = list( '\t'.join( list( entry.strip( ) for entry in [ line[ : 8 ], line[ 9 : 15 ], line[ 16 : 22 ], line[ 24 ], line[ 26 : 31 ], line[ 32 : 40 ], line[ 41 : 49 ], line[ 50 : ] ] ) ) for line in l_lines )
    df = pd.read_csv( StringIO( '\n'.join( l_csv_lines ) ), header = None, sep = '\t' )
    df.columns = [ 'UID', 'PID', 'PPID', 'C', 'STIME', 'TTY', 'TIME', 'CMD' ]
    return df


# In[ ]:


def OS_CPU_and_Memory_Usage( print_summary = True, return_summary_for_each_user = False, return_summary_for_each_program_of_each_user = False, ignore_processes_using_neglible_amount_of_resources = True ) :
    ''' Parse the output of a command line 'top -bn 1' and parse it into a DataFrame, and return the DataFrame '''
    l_line = os.popen( 'top -bn 1' ).read( ).split( '\n' )[ 6 : - 1 ]
    l_csv_line = list( '\t'.join( list( entry.strip( ) for entry in [ line[ : 6 ], line[ 7 : 15 ], line[ 16 : 19 ], line[ 20 : 23 ], line[ 24 : 31 ], line[ 32 : 38 ], line[ 39 : 45 ], line[ 46 ], line[ 48 : 53 ], line[ 54 : 58 ], line[ 59 : 68 ], line[ 69 : ] ] ) ) for line in l_line )
    df = pd.read_csv( StringIO( '\n'.join( l_csv_line ) ), sep = '\t' ).rename( columns = { '%CPU' : 'CPU', '%MEM' : 'MEM'  } ) # rename columns to make it compatible with pandas DataFrame
    if ignore_processes_using_neglible_amount_of_resources : df = df[ ( df.CPU > 0 ) | ( df.MEM > 0 ) ] # if 'ignore_processes_using_neglible_amount_of_resources' is true, remove processes using '0' CPU and '0' MEM.
    if print_summary : print( df[ [ 'CPU', 'MEM' ] ].sum( ) )
    if return_summary_for_each_user : return df[ [ 'USER', 'CPU', 'MEM' ] ].groupby( 'USER' ).sum( ).sort_values( 'MEM', ascending = False )
    elif return_summary_for_each_program_of_each_user : return df[ [ 'USER', 'COMMAND', 'CPU', 'MEM' ] ].groupby( [ 'USER', 'COMMAND' ] ).sum( ).sort_values( 'MEM', ascending = False )
    else : return df


# In[1]:


def OS_PIPELINE_Muitple_Running( arr_cmd_line, n_lines_at_a_time, dir_data = None, title = '', excute_cmd_line = False, wait_n_seconds_between_splited_jobs = 5, split = True, l_cmd_line_initialize = [ 'sleep 1', 'echo "Bookshelves PIPELINE Started"', 'date' ], l_server = None, server_replace_home_dir = False, shellscript_use_relative_path = False, dict_id_server_to_dir_home = { 'node210' : '/node210data/', "node01" : '/node01data/', "node200" : '/home/' } ) : # 2020-11-21 16:14:22 
    ''' # 2021-01-07 16:34:06 
    (Usage of 'split' argument is deprecated!) For a give list of cmd lines, execute 'n_lines_at_a_time' number of lines at one cycle, wait until all lines are completed, and start another cycle. 
    if cmd lines are not executed, write a shell script file at the directory given by dir_data if given.
    Additinally, report server time when each pipeline starts and when each command is completed. These information can be located by 'grep "Bookshelves"'
    if 'wait_n_seconds_between_splited_jobs' is given, wait n_seconds between initiating each job, so that resources can be evenly distributed and avoid bottle neck.
    Also, add shell command lines in 'l_cmd_line_initialize' at the start of each splitted file.
    
    l_server: should be a subset of [ 'node200', 'node01', 'node210' ]
              write multiple sets of shellscripts for running jobs in multiple servers
    server_replace_home_dir: if set to True, replace '/__current_server_home_dir__/' string with home directory of each server (useful for jobs requiring high file IO, such as RNA-Seq index or HMM database) 
    shellscript_use_relative_path : use relative path of child bash shell scripts when multiple shellscripts are run together. Should be set to True when running shellscripts in a docker environment.
    '''
    try : arr_cmd_line.dtype # preprocess the format of arr_cmd_line. 1) make sure arr_cmd_line is np.array datatype, and 2) dtype of the np.array is Object
    except : arr_cmd_line = np.array( arr_cmd_line, dtype = str ).astype( object )
    str_time_stamp = TIME_GET_timestamp( ) # get time stamp of current time to uniquely name the shellscript
    arr_cmd_line = arr_cmd_line.astype( object ) if arr_cmd_line.dtype is not np.dtype( 'O' ) else arr_cmd_line  # convert arr_cmd_line to object dtype array
    n_lines_at_a_time = min( n_lines_at_a_time, len( arr_cmd_line ) ) # set n_lines_at_a_time according to the given number of command lines
    if l_server is None: # if cmd_line is run in only one server # l_server = [ os.environ[ 'SERVER_NAME' ] ] # use current id_server by default
        arr_cmd_line_for_a_server = arr_cmd_line # retrieve cmd_line_for_each_server
        n_lines = len( arr_cmd_line_for_a_server )
        name_file = f'{title}__Muitple_Running_{len( arr_cmd_line_for_a_server )}_lines_{str_time_stamp}'
        for index_shell_script in range( n_lines_at_a_time ) :
            with open( f'{dir_data}{name_file}__splitted_{index_shell_script + 1}.sh', 'w' ) as file :
                arr_cmd_line_for_a_server_splitted = arr_cmd_line_for_a_server[ index_shell_script : : n_lines_at_a_time ]
                n_line_splitted = len( arr_cmd_line_for_a_server_splitted ) # add echo command line that report the progress of the given jobs for each splitted file 
                arr_percents = np.round( np.arange( 0, 100, 100 / n_line_splitted ) + 100 / n_line_splitted, 1 )
                arr_percents = arr_percents[ arr_percents <= 100 ] # remove occasional percent values above 100%
                arr_cmd_line_for_a_server_splitted = arr_cmd_line_for_a_server_splitted + '\necho "Bookshelves PIPELINE: ' + arr_percents.astype( str ).astype( object ) + ' % completed (line ' + np.arange( 1, n_line_splitted + 1 ).astype( str ).astype( object ) + ' at splitted file {})"\ndate\necho "Completed Command : '.format( index_shell_script + 1 ) + arr_cmd_line_for_a_server_splitted + '"' # echo the complete percentage of given jobs and echo the given command line name
                str_shell_script = '\n'.join( l_cmd_line_initialize ) + '\n' + '\n'.join( arr_cmd_line_for_a_server_splitted )
                if wait_n_seconds_between_splited_jobs : file.write( f'sleep {int( ( index_shell_script + 1 ) * wait_n_seconds_between_splited_jobs )}\n{str_shell_script}' ) # if 'wait_n_seconds_between_splited_jobs' is given, wait n_seconds between initiating each job, so that resources can be evenly distributed and avoid bottle neck.
                else : file.write( str_shell_script ) # add shell command lines in 'l_cmd_line_initialize' at the start of each splitted file.
        if n_lines_at_a_time > 1 :
            with open( f'{dir_data}{name_file}.sh', 'w' ) as file :
                if not shellscript_use_relative_path :
                    file.write( f"cd {dir_data}\n" )
                file.write( ' & '.join( list( f'bash {name_file if shellscript_use_relative_path else dir_data + name_file}__splitted_{index_shell_script + 1}.sh' for index_shell_script in range( n_lines_at_a_time ) ) ) + '\nwait' )
        else : os.rename( f'{dir_data}{name_file}__splitted_{index_shell_script + 1}.sh', dir_data + f'{name_file}.sh' )
    else : # if cmd_lines are run on multiple servers
        n_servers = len( l_server )
        for index_server, id_server in enumerate( l_server ) : # for each server
            arr_cmd_line_for_a_server = arr_cmd_line[ index_server : : n_servers ] # retrieve cmd_line_for_each_server
            dir_home_current_server = dict_id_server_to_dir_home[ id_server ]
            if server_replace_home_dir : # replace home directory in 'cmd_line' to the current server's home directory
                arr_cmd_line_for_a_server = np.array( list( cmd_line.replace( '/__current_server_home_dir__/', dir_home_current_server ) for cmd_line in arr_cmd_line_for_a_server ), dtype = object ) # modify the absolute path to match the dir_folder_data_server for each server
            n_lines = len( arr_cmd_line_for_a_server )
            n_lines_at_a_time_server = n_lines_at_a_time if n_lines >= n_lines_at_a_time else n_lines # if number of command lines is larger than current 'n_lines_at_a_time', set 'n_lines' as 'n_lines_at_a_time_server' for the current server
            name_file = f'{id_server}.{title}__Muitple_Running_{len( arr_cmd_line_for_a_server )}_lines_{str_time_stamp}'
            for index_shell_script in range( n_lines_at_a_time_server ) :
                with open( f'{dir_data}{name_file}__splitted_{index_shell_script + 1}.sh', 'w' ) as file :
                    arr_cmd_line_for_a_server_splitted = arr_cmd_line_for_a_server[ index_shell_script : : n_lines_at_a_time_server ]
                    n_line_splitted = len( arr_cmd_line_for_a_server_splitted ) # add echo command line that report the progress of the given jobs for each splitted file 
                    arr_percents = np.round( np.arange( 0, 100, 100 / n_line_splitted ) + 100 / n_line_splitted, 1 )
                    arr_percents = arr_percents[ arr_percents <= 100 ] # remove occasional percent values above 100%
                    arr_cmd_line_for_a_server_splitted = arr_cmd_line_for_a_server_splitted + '\necho "Bookshelves PIPELINE: ' + arr_percents.astype( str ).astype( object ) + ' % completed (line ' + np.arange( 1, n_line_splitted + 1 ).astype( str ).astype( object ) + ' at splitted file {})"\ndate\necho "Completed Command : '.format( index_shell_script + 1 ) + arr_cmd_line_for_a_server_splitted + '"' # echo the complete percentage of given jobs and echo the given command line name
                    str_shell_script = '\n'.join( l_cmd_line_initialize ) + '\n' + '\n'.join( arr_cmd_line_for_a_server_splitted )
                    if wait_n_seconds_between_splited_jobs : file.write( f'sleep {int( ( index_shell_script + 1 ) * wait_n_seconds_between_splited_jobs )}\n{str_shell_script}' ) # if 'wait_n_seconds_between_splited_jobs' is given, wait n_seconds between initiating each job, so that resources can be evenly distributed and avoid bottle neck.
                    else : file.write( str_shell_script ) # add shell command lines in 'l_cmd_line_initialize' at the start of each splitted file.
            if n_lines_at_a_time_server > 1 :
                with open( dir_data + '{name_file}.sh'.format( name_file = name_file ), 'w' ) as file :
                    if not shellscript_use_relative_path :
                        file.write( f"cd {dir_data}\n" )
                    file.write( ' & '.join( list( f'bash {name_file if shellscript_use_relative_path else dir_data + name_file}__splitted_{index_shell_script + 1}.sh' for index_shell_script in range( n_lines_at_a_time_server ) ) ) + '\nwait' )
            else : os.rename( f'{dir_data}{name_file}__splitted_1.sh', f'{dir_data}{name_file}.sh' )
    return f'{title}__Muitple_Running_{len( arr_cmd_line_for_a_server )}_lines_{str_time_stamp}' # return a filename to help to identify shellscripts written by the function


# In[ ]:


def OS_PIPELINE_Run_Until_Passed( l_popenargs, Func_check_passed, dict_args_Func_check_passed = dict( ) ) : # 2020-11-12 19:42:32 
    """ Run a process until its stdout and stderr outputs satisfy the given function 'Func_check_passed' with keyworded arguments 'dict_args_Func_check_passed'
     'Func_check_passed' : should return True if the output satisfy the conditions defining passed, and return False if the output does not satisfy the conditions. The function is called with the signature shown below:
                             Func_check_passed( stdout, stderr, ** dict_args_Func_check_passed )
     """
    str_id = UUID( ) # retrieve uuid of the current process

    dir_file_stdout = '/tmp/{uuid}.out.txt'.format( uuid = str_id ) # define stdout and stdin files
    dir_file_stderr = '/tmp/{uuid}.err.txt'.format( uuid = str_id )

    while True :
        with open( dir_file_stdout, 'w+' ) as fout : # excute and read std output and std errors of a process
            with open( dir_file_stderr, 'w+' ) as ferr :
                out = subprocess.call( l_popenargs, stdout = fout, stderr = ferr )
                fout.seek( 0 )
                stdout = fout.read( )
                ferr.seek( 0 ) 
                stderr = ferr.read( )
        if Func_check_passed( stdout, stderr, ** dict_args_Func_check_passed ) : break

    os.remove( dir_file_stdout ) # remove temporary files
    os.remove( dir_file_stderr )


# In[ ]:

def OS_Memory( ) :
    """
    # 2021-03-31 13:16:21 
    read and parse the memory of the system at '/proc/meminfo'
    """
    with open( '/proc/meminfo' ) as file :
        l_line = file.read( ).strip( ).split( '\n' )

    dict_mem = dict( )
    for line in l_line :
        l = line.split( )
        str_name_field = l[ 0 ][ : -1 ] # retrieve the name of the field # remove ':' character at the end
        dict_mem[ str_name_field ] = int( l[ 1 ] ) 
        if len( l ) > 2 : # if unit is available
            dict_mem[ f"{str_name_field}_unit" ] = l[ 2 ] # add unit to the dictionary
            
    return dict_mem

def OS_FILE_Combine_Files_in_order( l_dir_file, dir_newfile, overwrite_existing_file = False, delete_input_files = False, header = None, remove_n_lines = 0 ) : # 2020-07-20 11:47:29 
    ''' # 2021-01-04 14:51:40 
    combine contents of files in l_dir_file and write at dir_newfile. if header is given, append header (string type with \n at the end) at the front of the file. if 'remove_n_lines' > 0, remove n lines from each files.
    gzipped files are also supported. However, when input files and output files have mixed gzipped status, it will cause a TypeError '''
    if os.path.exists( dir_newfile ) and not overwrite_existing_file : print( "[OS_FILE_Combine_Files_in_order][ERROR] output file already exists" )
    elif len( l_dir_file ) == 0 : print( "[OS_FILE_Combine_Files_in_order][ERROR] given list of files is empty" )
    else : # if at least one input file is given (l_dir_file) and given directory of a newfile already exist (or overwriting is allowed)
        bool_flag_gzipped_output = dir_newfile[ - 3 : ] == '.gz' # set boolean flag for gzipped output file
        newfile = gzip.open( dir_newfile, 'wb' ) if bool_flag_gzipped_output else open( dir_newfile, 'w' ) # open a file that will contained combined contents
        if header : 
            header = header.decode( ) if not isinstance( header, ( str ) ) else header # convert header byte string to string if header is not a string type
            newfile.write( ( header.encode( ) if bool_flag_gzipped_output else header ) ) # write a header line to the output file 
        for dir_file in l_dir_file :
            bool_flag_gzipped_input = dir_file[ - 3 : ] == '.gz' # set boolean flag for gzipped input file
            with ( gzip.open( dir_file, 'rb' ) if bool_flag_gzipped_input else open( dir_file, 'r' ) ) as file :
                if remove_n_lines : 
                    for index in range( remove_n_lines ) : file.readline( )
                shutil.copyfileobj( file, newfile )
        newfile.close( )
        if delete_input_files :
            for dir_file in l_dir_file : os.remove( dir_file )


# ## Function Specific to a Project

# ### Function for TE_eQTL Project

# In[ ]:

def OS_Run( l_args, dir_file_stdout = None, dir_file_stderr = None, return_output = True, remove_default_output_files = True, stdout_binary = False ) :
    """ # 2021-03-30 19:41:16 
    Run a process and save stdout and stderr as a file.
    
    'return_output' : return the output as dictionary of strings
    'remove_default_output_files' : remove default stdout and stderr files containing the output of the process when 'dir_file_stdout' and 'dir_file_stderr' were not given.
    'stdout_binary' : set this flag to True if stdout is binary.
    """
    uuid_process = UUID( ) # set uuid of the process
    # define default stdout and stdin files and set approproate flags
    flag_dir_file_stdout_was_given = dir_file_stdout is not None
    flag_dir_file_stderr_was_given = dir_file_stderr is not None
    
    if not flag_dir_file_stdout_was_given :
        dir_file_stdout = f'/tmp/{uuid_process}.out.txt'
    if not flag_dir_file_stderr_was_given :
        dir_file_stderr = f'/tmp/{uuid_process}.err.txt'
    
    with open( dir_file_stdout, 'w+b' if stdout_binary else 'w+' ) as fout : # excute and read std output and std errors of a process
        with open( dir_file_stderr, 'w+' ) as ferr :
            out = subprocess.call( l_args, stdout = fout, stderr = ferr )
            fout.seek( 0 )
            stdout = fout.read( )
            ferr.seek( 0 ) 
            stderr = ferr.read( )
    # remove default output files
    if not flag_dir_file_stdout_was_given :
        os.remove( dir_file_stdout )
    if not flag_dir_file_stderr_was_given :
        os.remove( dir_file_stderr )
            
    return { 'stdout' : stdout, 'stderr' : stderr } if return_output else None


def TE_eQTL_Swarm_Plot( df_expr, df_genotype, index_entry, index_genotype ) :
    ''' Draw a swamplots showing differeing expressions in different genotypes '''
    df_eqtl_query = pd.DataFrame( dict( expr = df_expr.loc[ index_entry ], n_alt_alleles = df_genotype.loc[ int( index_genotype.rsplit( '_', 1 )[ 1 ] ) ] ) )    
    sns.swarmplot( data = df_eqtl_query, x = 'n_alt_alleles', y = 'expr' )
    MATPLOTLIB_basic_configuration( title = 'eQTL SwamPlot', x_label = 'number of alternative alleles of {}'.format( index_genotype ), y_label = 'detected expression of {}'.format( index_entry ) )


# ### Parsing function for file fomats frequently used in bioinformatics

# In[ ]:


def EMBL_Load_into_DataFrame( dir_file ) :
    ''' Read EMBL format annotation of a given file (should be unzipped) and parse annotations into a DataFrame '''
    l_dict_record = list( )
    with open( dir_file, 'r' ) as file :
        while True :
            dict_record = dict( )
            while True :
                line = file.readline( )
                if line == '//\n' or len( line ) == 0 : break
                str_tag, str_content = line[ : 2 ], line[ 5 : - 1 ].strip( )
                if str_tag != 'XX' and len( str_content ) > 0 :
                    character_for_joining = ' '
                    if str_tag == '  ' : 
                        str_tag, character_for_joining = 'Sequence', ''
                        str_content = STR.Replace_a_character_every_n_characters( str_content[ : 65 ], n_characters = 11, replacement_characters = '' ).strip( )
                    if str_content[ - 1 ] == '.' or str_content[ - 1 ] == ';' : str_content = str_content[ : - 1 ] # if content ends with ';' or '.', remove the characters.
                    if str_tag in dict_record : dict_record[ str_tag ] = dict_record[ str_tag ] + character_for_joining + str_content # combine content with same str_tag together
                    else : dict_record[ str_tag ] = str_content
            if len( line ) == 0 : break
            l_dict_record.append( dict_record )
    return pd.DataFrame( l_dict_record )


# ### Function for NGS data

# In[ ]:


def SEQ_Trim_Invalid_Residue( seq, invalid_residues = 'JOUXZB' ) : # 2020-07-14 16:25:34 
    """ Trim a given set of invalid residues (in a string format) from both ends from given sequence. If a stretch of invalid residues exist inside a trimmed sequence, retrieve either 3' valid sequence or 5' valid sequence depending on the length of the valid sequences. """
    l_invalid_residues = list( invalid_residues )
    length_seq = len( seq )
    char_representing_all_invalid_residue = l_invalid_residues[ 0 ]
    for char_invalid_residue in l_invalid_residues[ 1 : ] : seq = seq.replace( char_invalid_residue, char_representing_all_invalid_residue ) # replace all invalid_residues to a character representing all invalid residues
    n_invalid_residue_at_5_prime, n_invalid_residue_at_3_prime = 0, 0 # initialize number of N bases at either ends
    while n_invalid_residue_at_5_prime < length_seq and seq[ n_invalid_residue_at_5_prime ] == char_representing_all_invalid_residue : n_invalid_residue_at_5_prime += 1 # find out where N base ends at either end
    if n_invalid_residue_at_5_prime == length_seq : return '' # if a sequence contains only char_representing_all_invalid_residue base, return an empty fastq record
    while seq[ length_seq - n_invalid_residue_at_3_prime - 1 ] == char_representing_all_invalid_residue : n_invalid_residue_at_3_prime += 1
    seq_trimmed = seq[ n_invalid_residue_at_5_prime : length_seq - n_invalid_residue_at_3_prime ] # trim a stretch of N bases at either end
    if char_representing_all_invalid_residue in seq_trimmed : # if char_representing_all_invalid_residue base exist inside a trimmed sequence, retrieve either 3' non-N base stretch or 5' non-N base stretch depending on the length of the non-N base stretchs.
        n_non_invalid_residue_at_5_prime = seq_trimmed.find( char_representing_all_invalid_residue )
        n_non_invalid_residue_at_3_prime = seq_trimmed[ : : - 1 ].find( char_representing_all_invalid_residue )
        if n_non_invalid_residue_at_5_prime >= n_non_invalid_residue_at_3_prime : seq_trimmed = seq_trimmed[ : n_non_invalid_residue_at_5_prime ]
        else : seq_trimmed = seq_trimmed[ - n_non_invalid_residue_at_3_prime : ]
    return seq_trimmed


# In[ ]:


def SAM_Header_to_Chrom_Length_DataFrame( l_header ) :
    """  prepare dataframe for calculating global coordinate by using headers from SAM output file  """ 
    df = pd.DataFrame( np.array( list( list( entry[ 3 : ] for entry in sorted( header[ 4 : ].split( '\t' ) )[ : : -1 ] ) for header in l_header if '@SQ' in header ) ), columns = [ 'seqname', 'seqlength' ] )
    df.seqlength = df.seqlength.astype( int )
    position_global = 0 
    l_position_global = list( )
    for seqlength in df.seqlength :
        l_position_global.append( position_global )
        position_global += seqlength
    df[ 'start_position_global' ] = l_position_global
    return df


# In[ ]:


def NGS_CIGAR_Iterate_a_CIGAR_String( str_cigar ) :
    ''' Iterate through a CIGAR String. yield an integer and operator, which is one of MIDNSHP=X '''
    str_number = ''
    for char in str_cigar :
        if char in 'MIDNSHP=X' :
            yield int( str_number ), char
            str_number = ''
        else :
            str_number += char


# In[ ]:


def NGS_CIGAR_Retrive_List_of_Mapped_Segments( cigar, pos_global ) :
    ''' return l_seq and int_total_aligned_length for given cigar string and pos_global '''
    l_seg, start, int_aligned_length, int_total_aligned_length = list( ), pos_global, 0, 0
    for length, operation in NGS_CIGAR_Iterate_a_CIGAR_String( cigar ) :
        if operation in 'MD=X' : int_aligned_length += length
        elif operation == 'N' : # if splice junction appears, split the region and make a record
            l_seg.append( ( start, start + int_aligned_length - 1 ) ) # set the end position
            start = start + int_aligned_length + length # set the next start position
            int_total_aligned_length += int_aligned_length # increase total aligned length
            int_aligned_length = 0
    if int_aligned_length > 0 : 
        l_seg.append( ( start, start + int_aligned_length - 1 ) )
        int_total_aligned_length += int_aligned_length
    return l_seg, int_total_aligned_length


# In[ ]:


def NGS_GTF_Parse_Attribute( record ) :
    ''' Parse GTF attribute into a dictionary '''
    return ast.literal_eval( '{' + ', '.join( list( '"' + entry.replace( ' "', '" : "' ) for entry in record.strip( )[ : -1 ].split( '; ' ) ) ) + '}' )


# In[ ]:


def NGS_SEQ_Generate_Kmer( seq, window_size ) :
    """ 
    # 2021-02-20 15:14:13 
    generate a list of Kmers from the sequence with the given window_size  """
    return list( seq[ i : i + window_size ] for i in range( 0, len( seq ) - window_size + 1, 1 ) )


# In[ ]:


def NGS_SEQ_Reverse_Complement( seq ) :
    ''' # 2021-02-04 11:47:19 
    Return reverse complement of 'seq' '''
    dict_dna_complement = { "A" : 'T', "T" : 'A', "C" : 'G', "G" : 'C', "N" : 'N', "-" : '-' }
    return ''.join( list( dict_dna_complement[ base ] for base in seq ) )[ : : -1 ]


# In[ ]:


# 2020-05-29 10:53:13 
dict_NGS__aa_to_codon = { 'I' : [ 'ATT', 'ATC', 'ATA' ], 'L' : [ 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG' ], 'V' : [ 'GTT', 'GTC', 'GTA', 'GTG' ], 'F' : [ 'TTT', 'TTC' ], 'M' : [ 'ATG' ], 'C' : [ 'TGT', 'TGC' ], 'A' : [ 'GCT', 'GCC', 'GCA', 'GCG' ], 'G' : [ 'GGT', 'GGC', 'GGA', 'GGG' ], 'P' : [ 'CCT', 'CCC', 'CCA', 'CCG' ], 'T' : [ 'ACT', 'ACC', 'ACA', 'ACG' ], 'S' : [ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ], 'Y' : [ 'TAT', 'TAC' ], 'W' : [ 'TGG' ], 'Q' : [ 'CAA', 'CAG' ], 'N' : [ 'AAT', 'AAC' ], 'H' : [ 'CAT', 'CAC' ], 'E' : [ 'GAA', 'GAG' ], 'D' : [ 'GAT', 'GAC' ], 'K' : [ 'AAA', 'AAG' ], 'R' : [ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ], '*' : [ 'TAA', 'TAG', 'TGA' ] }
dict_NGS__codon_to_aa = dict( )
for aa in dict_NGS__aa_to_codon :
    l_codon = dict_NGS__aa_to_codon[ aa ]
    for codon in l_codon : dict_NGS__codon_to_aa[ codon ] = aa
        
def NGS_SEQ_Translate( seq, start = 0, append_stop_codon_to_returned_seq = False, print_message = True, return_error_value_if_stop_codon_not_found = True ) :
    ''' Translate nucleotide sequence until the given seq is end or stop codons are encountered.
    start = 0 : start position for potential protein sequence (for example, start of ATG start codon.) '''
    if type( seq ) is float : return np.nan # if input is invalid, return np.nan
    length_seq = len( seq ) # it seems the time complexity of the len(str) method is O(1)
    l_aa = list( )
    while True :
        if start + 3 > length_seq : 
            if print_message : print( 'Stop Codon Not Found' )
            if return_error_value_if_stop_codon_not_found : return np.nan
            else : break
        try : aa = dict_NGS__codon_to_aa[ seq[ start : start + 3 ] ] # try to translat the sequence
        except : 
            if return_error_value_if_stop_codon_not_found : return np.nan # if an error occur during translating and stop codon was not found, return np.nan if 'return_error_value_if_stop_codon_not_found'
            else : break
        if aa == '*' :
            if append_stop_codon_to_returned_seq : l_aa.append( aa )
            break
        else : l_aa.append( aa )
        start += 3
    return ''.join( l_aa )


# In[ ]:


# 2020-05-29 14:28:22 
def NGS_SEQ_ORF_Find_All_Methionine_ORFs( seq, return_list = False ) :
    ''' Find all ORFs of a given sequence, and return it as a dataframe ( columns = [ 'start', 'end', 'translated_protein_seq' ], use 0-based coordinate for start and end ). 
    return found ORFs as 2D list if 'return_list' is set to True. '''
    l_l_value = [ [ 0, 0, np.nan ] ] # 2D array with a dummy starting ORF
    for int_atg_start_pos in STR.Find_all( seq, 'ATG' ) : # find all positions of subsequences starting with ATG (methionine start codon)
        prot_seq = NGS_SEQ_Translate( seq, int_atg_start_pos, append_stop_codon_to_returned_seq = True, return_error_value_if_stop_codon_not_found = False, print_message = False )
        end_pos = int_atg_start_pos + 3 * len( prot_seq ) # retrive end position of CDS, including a stop codon if detected.
        if l_l_value[ - 1 ][ 1 ] == end_pos : continue # if current ORF shared the same stop position (and with downstream start position), do not include the ORF in the output, since current ORF is a fragment of the previous ORF
        else : l_l_value.append( [ int_atg_start_pos, end_pos, prot_seq ] ) # use 0-based coordinate
    return l_l_value[ 1 : ] if return_list else pd.DataFrame( l_l_value[ 1 : ], columns = [ 'start', 'end', 'translated_protein_seq' ] ) # return found ORFs without the dummy ORF at the start of the 2D array.
# 2020-05-29 14:28:22 
def NGS_SEQ_ORF_Find_All_Methionine_ORFs_on_Both_Strands( seq, use_1_based_coordinate = True ) :
    ''' Find all ORFs start with methionine on both strand of the given sequence by using 'NGS_SEQ_ORF_Find_All_Methionine_ORFs'. 
    use 1-based-coordinate for start and end positions if 'use_1_based_coordinate' is set to True '''
    df_plus = NGS_SEQ_ORF_Find_All_Methionine_ORFs( seq ) # find all ORFs in the sequence and the reverse complement of the sequence
    df_minus = NGS_SEQ_ORF_Find_All_Methionine_ORFs( NGS_SEQ_Reverse_Complement( seq ) )
    df_plus[ 'strand' ] = '+'
    df_minus[ 'strand' ] = '-'
    int_length_seq = len( seq ) # convert coordinate of reverse complement of a sequence into that of the sequence
    s_start = int_length_seq - df_minus.end
    s_end = int_length_seq - df_minus.start
    df_minus[ 'start' ] = s_start
    df_minus[ 'end' ] = s_end

    df_orf = pd.concat( [ df_plus, df_minus ], ignore_index = True ) # combine ORFs from plus strand and the minus strand
    df_orf.sort_values( 'start', ignore_index = True, inplace = True )
    if use_1_based_coordinate : df_orf.start += 1 # use 1-based-coordinate for start and end positions if 'use_1_based_coordinate' is set to True
    return df_orf


# In[ ]:


def NGS_SEQ_Trim_PolyA( seq, from_3_prime_end = True, return_length_of_polyA = False, str_repeated_base = 'A', int_lookup_window = 3 ) :
    ''' Trim PolyA sequence from the 3' end of a given sequence (DNA sequence in upper characters) if 'from_3_prime_end' is True or trim polyA from 5' end of the sequence if 'from_3_prime_end' is False
    if the program encounter base other than 'A', lookup next 'int_lookup_window' number of bases and check weather they are consecutive 'A' bases. '''
    if len( seq ) == 0 : # if the given seq is empty 
        if return_length_of_polyA : return 0
        else : return ''
    str_repeated_bases_in_lookup_window = str_repeated_base * int_lookup_window
    seq_polyA_at_5_prime_end = seq[ : : -1 ] if from_3_prime_end else seq
    for index, str_base in enumerate( seq_polyA_at_5_prime_end ) :
        if str_base != str_repeated_base :
            if seq_polyA_at_5_prime_end[ index + 1 : index + 1 + int_lookup_window ] == str_repeated_bases_in_lookup_window : continue # even though sequence is not 'A', look up the next bases, and if the next consequtive three bases are 'A', consider the character as heterogenity of PolyA and countinue counting number of 'A's in PolyA 
            else : break
    if return_length_of_polyA : return index
    else : return ( seq[  : - index ] if index > 0 else seq ) if from_3_prime_end else seq[ index : ]


# In[ ]:


dict_NGS_encoding_seq_to_int = { '' : 0,  'A' : 1 , 'C' : 2 , 'G' : 3 , 'T' : 4  }
def NGS_SEQ_Encode_to_integer( seq, reverse_order_during_encoding = True ) :
    ''' convert sequence string ( should be in upper case ) with maximum length 27bp into integer (64 bit integer) in base-5 numeral system ('A' = 1, 'C' = 2, 'G' = 3, 'T' = 4). For algorithmical reason, encoded sequence is in an inverse order. '''
    int_encoded_seq = 0
    if not reverse_order_during_encoding : seq = seq[ : : -1 ] 
    for index, base in enumerate( seq ) :
        int_encoded_seq += dict_NGS_encoding_seq_to_int[ base ] * 5 ** index
    return int_encoded_seq

dict_NGS_decoding_int_to_seq = { 0 : '', 1 : 'A', 2 : 'C', 3 : 'G', 4 : 'T' }

def NGS_SEQ_Decode_from_integer( int_encoded_seq, reverse_order_during_encoding = False ) :
    ''' convert integer (64 bit integer) in base-5 numeral system ('A' = 1, 'C' = 2, 'G' = 3, 'T' = 4) into sequence string ( in upper case ). For algorithmical reason, decoded sequence is in an inverse order. '''
    l_decoded_seq = list( )
    while True :
        decoded_base = dict_NGS_decoding_int_to_seq[ int_encoded_seq % 5 ]
        if not decoded_base : break
        l_decoded_seq.append( decoded_base )
        int_encoded_seq = int( int_encoded_seq / 5 )
    seq = ''.join( l_decoded_seq )
    return seq if reverse_order_during_encoding else seq[ : : -1 ]


# In[ ]:


def NGS_SEQ_Calculate_Simple_Repeat_Proportion_in_a_read( seq, int_len_kmer = 4, int_kmer_count_threshold_for_repeat = 5 ) :
    ''' Calculate proportion of simple repeat in a read. simple repeat is detected by counting kmer of the given length and retrive kmer with count higher than the given threshold.  '''
    len_seq = len( seq )
    n_kmers = len_seq - int_len_kmer + 1
    dict_kmer_count = COUNTER( seq[ index : index + int_len_kmer ] for index in range( n_kmers ) ) # count kmers of the given length
    int_kmer_count_repeat = 0 # count kmers that have count values above the given threshold
    for kmer in dict_kmer_count :
        int_kmer_count = dict_kmer_count[ kmer ]
        if int_kmer_count > int_kmer_count_threshold_for_repeat : int_kmer_count_repeat += int_kmer_count
    float_ratio_of_repeat = int_kmer_count_repeat / n_kmers
    return float_ratio_of_repeat


# In[ ]:


def FASTA_Read( dir_fasta, print_message = False, remove_space_in_header = False, return_dataframe = False, parse_uniprot_header = False ) : # 2020-08-21 14:38:09 
    ''' for a file-like object of file of 'dir_fasta' directory, parse the content into a dictionary. 'remove_space_in_header' option remove space in the header line (required for ABySS output)
    'parse_uniprot_header': if set to True and 'return_dataframe' is set to True, parse uniprot sequence header into [ 'accession', 'description', 'uniprot_source', 'uniprot_acc', 'uniprot_name' ] '''
    dict_header_to_seq = dict( )
    dict_duplicated_header_count = dict( )
    bool_flag_input_gzipped = False
    if hasattr( dir_fasta, 'readline' ) : file = dir_fasta # if file-like object was given instead of dir_fasta, use dir_fasta as file
    else : 
        bool_flag_input_gzipped = dir_fasta[ - 3 : ] == '.gz' 
        file = gzip.open( dir_fasta, 'rb' ) if bool_flag_input_gzipped else open( dir_fasta, 'r' ) # open file of dir_fasta depending on the detected gzipped status
    line = ( file.readline( ).decode( ) if bool_flag_input_gzipped else file.readline( ) )[ : -1 ]
    while True :
        str_header = line
        l_seq = list( )
        while True :
            line = ( file.readline( ).decode( ) if bool_flag_input_gzipped else file.readline( ) )[ : -1 ]
            if not line or line[ : 1 ] == '>' : break
            l_seq.append( line )
        header, seq = str_header[ 1 : ], ''.join( l_seq )
        if header in dict_duplicated_header_count : dict_duplicated_header_count[ header ] += 1 # handle sequences with duplicated headers
        else : dict_duplicated_header_count[ header ] = 1
        if dict_duplicated_header_count[ header ] > 1 in dict_header_to_seq : dict_header_to_seq[ header + '_dup_{index}'.format( index = dict_duplicated_header_count[ header ] ) ] = ''.join( l_seq ) # if current fastq header in already exists, add '_dup_{index}' to the header to make it unique 
        else : dict_header_to_seq[ header ] = ''.join( l_seq )
        if not line : break
    file.close( ) # close the file
    if len( dict_header_to_seq ) == 1 and tuple( dict_header_to_seq ) == tuple( [ '' ] ) : dict_header_to_seq = dict( ) # if 'dict_header_to_seq' is empty, return an empty dictionary as 'dict_header_to_seq'
    if print_message : print( pd.Series( dict_duplicated_header_count ) )
    if remove_space_in_header : dict_header_to_seq = dict( ( header.replace( ' ', '_' ), dict_header_to_seq[ header ] ) for header in dict_header_to_seq ) #  remove space in the header line
    if return_dataframe : # return parsed fasta file as a dataframe 
        df = pd.Series( dict_header_to_seq ).reset_index( )
        df.columns = [ 'header', 'seq' ]
        df[ 'length' ] = df.seq.apply( len )
        if parse_uniprot_header : # if set to True and 'return_dataframe' is set to True, parse uniprot sequence header into [ 'accession', 'description', 'uniprot_source', 'uniprot_acc', 'uniprot_name' ]
            l_l_value = list( )
            for header in df.header.values :
                accession, description = header.split( ' ', 1 )
                uniprot_source, uniprot_acc, uniprot_name = accession.split( '|' )
                l_l_value.append( [ accession, description, uniprot_source, uniprot_acc, uniprot_name ] )
            df = df.join( pd.DataFrame( l_l_value, columns = [ 'accession', 'description', 'uniprot_source', 'uniprot_acc', 'uniprot_name' ] ) )
        return df
    else : return dict_header_to_seq


# In[ ]:


def FASTA_Write( dir_fasta, dict_fasta = None, l_id = None, l_seq = None, overwrite_existing_file = False ) :
    ''' write fasta file at the given directory with dict_fastq (key = fasta_header, value = seq) or given list of id (fasta_header) and seq '''
    if os.path.exists( dir_fasta ) and not overwrite_existing_file : 
        print( 'the file already exists' )
        return -1
    if dict_fasta is not None : # if dict_fasta was given.
        with open( dir_fasta, 'w' ) as file : file.write( ''.join( list( ">{}\n{}\n".format( name, STR.Insert_characters_every_n_characters( dict_fasta[ name ], 60, insert_characters = '\n' ) ) for name in dict_fasta ) ) )
    else : # if l_id and l_seq were given.
        if type( l_id ) is str and type( l_seq ) is str : # if only one sequence and sequence name were given
            with open( dir_fasta, 'w' ) as file : file.write( ">{}\n{}\n".format( l_id, STR.Insert_characters_every_n_characters( l_seq, 60, insert_characters = '\n' ) ) )
        else : # if list of sequences and sequence names were given
            with open( dir_fasta, 'w' ) as file : file.write( ''.join( list( ">{}\n{}\n".format( name, STR.Insert_characters_every_n_characters( seq, 60, insert_characters = '\n' ) ) for name, seq in zip( l_id, l_seq ) ) ) )


# In[ ]:


def FASTA_DNA_break_scaffolds( dict_fasta, remove_space_in_header = True ) :
    ''' preprocess assembled contig before flye subassembly. replace space characers in the fasta header with '_' and break scaffolds containing 'NNN' (unknown length of gaps) into contigs '''
    dict_fasta = FASTA_Read( dict_fasta ) if isinstance( dict_fasta, ( str ) ) else dict_fasta # read fasta file if path (string) is given, or use given dictionary as a dict_fasta
    if remove_space_in_header : dict_fasta = dict( ( header.replace( ' ', '_' ), dict_fasta[ header ] ) for header in dict_fasta ) # if 'remove_space_in_header' was set to True, replace space characers in the fasta header with '_'
    dict_fasta_scaffolds_broken = dict( )
    for header in dict_fasta :
        seq = dict_fasta[ header ].upper( ) # convert lower character 'n' with 'N'
        if 'N' in seq : # if 'N' is contained in the sequence, break scaffolds into contigs and modify fasta header 
            for index, seq_broken_scaffolds in enumerate( seq.replace( 'N', ' ' ).split( ) ) : dict_fasta_scaffolds_broken[ header + '_broken_scaffold_{index}'.format( index = index ) ] = seq_broken_scaffolds
        else : dict_fasta_scaffolds_broken[ header ] = seq
    return dict_fasta_scaffolds_broken


# In[ ]:


def FASTA_Filter( dict_fasta, int_thres_length ) :
    ''' filter out fasta records smaller than the given threshold, 'int_thres_length'  '''
    dict_fasta = FASTA_Read( dict_fasta ) if isinstance( dict_fasta, ( str ) ) else dict_fasta # read fasta file if path (string) is given, or use given dictionary as a 'dict_fasta'
    for header in set( header for header in dict_fasta if int_thres_length > len( dict_fasta[ header ] ) ) : dict_fasta.pop( header ) # retrive set of headers with sequence length smaller than 'int_thres_length' and remove them from the dictionary
    return dict_fasta


# In[ ]:


# 2020-05-29 23:09:52 
def FASTA_Assembly_Stats( dict_fasta ) :
    ''' return assembly stats for the given fasta file (fasta file containing assembed contigs)  '''
    dict_fasta = FASTA_Read( dict_fasta ) if isinstance( dict_fasta, ( str ) ) else dict_fasta # read fasta file if path (string) is given, or use given dictionary as a 'dict_fasta'
    int_n_records = len( dict_fasta )
    arr_length = np.zeros( int_n_records )
    for index, header in enumerate( dict_fasta ) :
        arr_length[ index ] = len( dict_fasta[ header ] )
    arr_length.sort( )
    arr_length = arr_length[ : : -1 ] # reverse the order so that largest contigs are at the front
    arr_length_percentage = arr_length / arr_length.sum( ) * 100    
    arr_length_percentage_cumulative = np.zeros( int_n_records )
    arr_length_percentage_cumulative[ 0 ] = arr_length_percentage[ 0 ]
    for index, percentage in enumerate( arr_length_percentage[ 1 : ], 1 ) :
        arr_length_percentage_cumulative[ index ] = arr_length_percentage_cumulative[ index - 1 ] + percentage # calculate cumulative length
    l_values = [ arr_length.sum( ), arr_length[ 0 ] ] # retrive total_length and length of the longest contig
    for int_thres_percentage in np.arange( 10, 100, 10 ) :
        mask = arr_length_percentage_cumulative < int_thres_percentage # mask for contigs below thres_percentage when cumulated percentage is calculated by using larger contigs first.
        l_values.extend( [ arr_length[ np.where( np.diff( mask ) )[ 0 ][ 0 ] ], mask.sum( ) ] ) # calculate assembly-stats
    l_values.extend( [ arr_length.min( ), len( mask ) ] ) # append values for N100
    s_assembly_stats = pd.Series( l_values, index = [ 'total_length', 'longest', 'N10', 'N10n', 'N20', 'N20n', 'N30', 'N30n', 'N40', 'N40n', 'N50', 'N50n', 'N60', 'N60n', 'N70', 'N70n', 'N80', 'N80n', 'N90', 'N90n', 'N100', 'N100n' ] ).astype( int )
    return s_assembly_stats


# In[ ]:


def GENOME_Add_Global_Position( df, col_chrom, col_start, col_end ) :
    ''' Add global position to GTF-like annotation dataframe '''
    if 'dict_chr_to_chr_start_global' not in globals( ) : globals( )[ 'dict_chr_to_chr_start_global' ] = DB_Load_a_data_object( 'df_anno_GRCh38_Ensembl__chrom_length', dataset = 'TE_eQTL' ).set_index( 'seqname' ).start_position_global.to_dict( ) # load 'dict_chr_to_chr_start_global' if it has not been loaded.
    MAP.dict_a2b = dict_chr_to_chr_start_global
    col_chrom_start_genome = col_chrom + '_start_genome'
    df[ col_chrom_start_genome ] = df[ col_chrom ].apply( MAP.Map_a2b )
    df[ col_start + '_genome' ] = df[ col_start ] + df[ col_chrom_start_genome ]
    df[ col_end + '_genome' ] = df[ col_end ] + df[ col_chrom_start_genome ]
    return df.drop( columns = [ col_chrom_start_genome ] )


# In[ ]:


def MSA_Calculate_Base_Frequency( l_seq, show_heatmap = False, cmap = 'viridis' ) :
    ''' return dataframe containing base frequencies for each position in the list of aligned sequences. Either RNA/DNA/Protein Sequence can be given as an input '''
    arr_seq = np.array( list( list( str_aligned_seq ) for str_aligned_seq in l_seq ), dtype = object ) # build an array of bases 
    df_base_frequency = pd.DataFrame( list( LIST_COUNT( arr_base_at_position, duplicate_filter = None ) for arr_base_at_position in arr_seq.T ) ).T.fillna( 0 ) / len( arr_seq ) * 100
    if show_heatmap : # show heatmap of base frequency if 'show_heatmap' is True
        plt.imshow( df_base_frequency.values, cmap = cmap, aspect = 7 )
        plt.yticks( np.arange( len( df_base_frequency ) ), labels = df_base_frequency.index.values )
        MATPLOTLIB_basic_configuration( show_grid = False, y_lim = dict( top = - 0.5, bottom = len( df_base_frequency ) - 0.5 ), show_colorbar = True )
    return df_base_frequency


# ## Functions for working with NEWICK tree (lineage tree)

# In[ ]:


def NEWICK_Parse_by_Recursion( str_newick_tree, default_branch_length = 1, verbose = False ) :
    ''' # 2020-12-10 19:51:25 
    parse newick tree by recursion
    'default_branch_length' : default branch length when only node name was given (e.g. '(A,B,(C,D));' )
     '''
    count_layer_parentheses = 0
    pos = 0
    len_str_newick_tree = len( str_newick_tree )
    str_node = ''
    dict_tree = dict( )
    while pos < len_str_newick_tree :
        character = str_newick_tree[ pos ]
        if count_layer_parentheses == 0 and ( character == ',' or pos == len_str_newick_tree - 1 ) : # parse a node 
            if pos == len_str_newick_tree - 1 :
                if character != ';' :
                    str_node += character
                else :
                    if verbose :
                        print( "the end of the newick tree" )
            str_node = str_node.strip( ) # remove empty spaces
            if len( str_node ) > 0 and str_node[ 0 ] == '(' : # if node not empty and is a nested newick tree
                str_nested_newick_tree, str_name_and_branch_length_of_node = str_node.rsplit( ')', 1 )
                str_nested_newick_tree = str_nested_newick_tree[ 1 : ] # remove the initial '(', too, to remove it from the nested structure.
                if verbose : 
                    print( "nested tree is " + str_nested_newick_tree[ : 10  ] )
                dict_subtree = NEWICK_Parse_by_Recursion( str_nested_newick_tree, default_branch_length = default_branch_length, verbose = verbose )
            else : # if node is a single node
                str_name_and_branch_length_of_node = str_node
                dict_subtree = dict( )
                
            if ':' in str_name_and_branch_length_of_node :
                str_name_node, str_float_branch_length = str_name_and_branch_length_of_node.split( ':' )
                float_branch_length = default_branch_length if len( str_float_branch_length ) == 0 else float( str_float_branch_length ) # assign default branch length if the entry is empty
            else :
                str_name_node = str_name_and_branch_length_of_node
                float_branch_length = default_branch_length
            if len( str_name_node ) == 0 :
                str_name_node = UUID( ) # assign random node name if no node name was given.
            if str_name_node in dict_tree :
                str_name_node += '_' + UUID( ) # if current node name is not unique, add unique identifier (random barcode) to the node name
            dict_tree[ str_name_node ] = dict( name_node = str_name_node, branch_length = float_branch_length, dict_subtree = dict_subtree ) # add node to the tree
            str_node = '' # empty the string containing a node 
        else : # add a character to the string representing a node
            str_node += character 
        # adjust the level lf layers of parentheses
        if character == '(' :
            count_layer_parentheses += 1
        elif character == ')' :
            count_layer_parentheses -= 1
            
        pos += 1
    return dict_tree

def NEWICK_Traverse__Get_All_Node_Names( dict_newick_tree ) :
    ''' # 2020-12-10 19:51:11 
    get all node_names for the given parsed dict_newick_tree '''
    l_name_node = list( )
    for name_node in dict_newick_tree :
        l_name_node.append( name_node ) # add name_node to the list
        record = dict_newick_tree[ name_node ]
        dict_subtree = record[ 'dict_subtree' ]
        l_name_node.extend( NEWICK_Traverse__Get_All_Node_Names( dict_subtree ) )
    return l_name_node

def NEWICK_Get_Coordinates_of_Nodes_and_Edges( dict_newick_tree, float_space_between_neighbors = 0.01 ) :
    ''' # 2020-12-10 20:56:09 
     ===== Returns =====
    return 'arr_coordinates_subtree' and 'arr_edges_subtree'
    'arr_coordinates_subtree' : shape = ( n_nodes, 3 ). Each row represent x, y coordinates, and name_node, where a root is located at the left (x = 0) and the leaves at the right side.
    'arr_edges_subtree' : shape ( n_edges, 2 ). Contains list of edges, where each row (edge) represent 'node_name' of parent and child
    '''
    if isinstance( dict_newick_tree, ( str ) ) : # if a newick tree is given as a string, parse the string first 
        dict_newick_tree = NEWICK_Parse_by_Recursion( dict_newick_tree )
    float_positon_top_of_tree_previous = 0 # previous position of the top of all tree (y coordinates)
    float_positon_top_of_tree = 0 # current position of the top of all tree (y coordinates)
    l_arr_coordinates = list( ) # list of coordinates
    l_arr_edges = list( ) # list of edges
    for name_node in dict_newick_tree :
        record = dict_newick_tree[ name_node ]
        dict_subtree = record[ 'dict_subtree' ]
        if len( dict_subtree ) > 0 : # if current node contains a subtree
            arr_coordinates_subtree, arr_edges_subtree = NEWICK_Get_Coordinates_of_Nodes_and_Edges( dict_subtree )
            n_nodes = len( arr_coordinates_subtree ) # retrieve the number of nodes
            float_branch_length = record[ 'branch_length' ]
            arr_coordinates_subtree[ :, 0 ] += float_branch_length # push subtree to right by 'float_branch_length'
            arr_coordinates_subtree[ :, 1 ] += float_positon_top_of_tree + float_space_between_neighbors # lift subtree to the top with a given space between subtrees
            float_positon_top_of_tree_previous = float_positon_top_of_tree_previous # update 'float_positon_top_of_tree_previous'
            float_positon_top_of_tree = np.max( arr_coordinates_subtree[ :, 1 ] ) if n_nodes > 20 else max( arr_coordinates_subtree[ :, 1 ] ) # update the topest position of the tree (max is faster than np.max when number of nodes is small due to minimal overhead)
            l_arr_coordinates.append( arr_coordinates_subtree )
            if len( arr_edges_subtree ) > 0 : # if edges of subtree is not empty (there is another subtree inside the subtree), add the edges to the list of edges
                l_arr_edges.append( arr_edges_subtree )
            l_edges = list( [ name_node_child, name_node ] for name_node_child in dict_subtree )
            l_arr_edges.append( np.array( l_edges, dtype = object ) ) # add current node to the list of coordinates
        else : # if there is no subtree
            float_positon_top_of_tree_previous = float_positon_top_of_tree # update the topest position of the tree
            float_positon_top_of_tree += float_space_between_neighbors
        l_arr_coordinates.append( np.array( [ 0, ( float_positon_top_of_tree + float_positon_top_of_tree_previous ) / 2, name_node ], dtype = object ) ) # add current node to the list of coordinates
    return ( np.vstack( l_arr_coordinates ) if len( l_arr_coordinates ) > 1 else l_arr_coordinates ), ( np.vstack( l_arr_edges ) if len( l_arr_edges ) > 0 else l_arr_edges )


# ## Function for working with K-mers

# In[ ]:


def Kmer_Calculater_Kmer_Frequency( dict_fasta, int_length_kmer = 6 ) :
    ''' Calculate Kmer frequency for a dictionary containing sequences and return a dataframe containing Kmer frequency '''
    dict_seq_id_to_dict_kmer_count, dict_seq_id_to_length = dict( ), dict( )
    for seq_id in dict_fasta :    
        seq = dict_fasta[ seq_id ] 
        dict_kmer_count = dict( ) # count Kmer for each sequence
        for int_start_kmer in range( len( seq ) - int_length_kmer + 1 ) : 
            kmer = seq[ int_start_kmer : int_start_kmer + int_length_kmer ]
            if kmer in dict_kmer_count : dict_kmer_count[ kmer ] += 1
            else : dict_kmer_count[ kmer ] = 1
        dict_seq_id_to_dict_kmer_count[ seq_id ] = dict_kmer_count
        dict_seq_id_to_length[ seq_id ] = len( seq )
    df = pd.DataFrame( dict_seq_id_to_dict_kmer_count ).fillna( 0 )
    s = pd.Series( dict_seq_id_to_length ) - int_length_kmer + 1 # calculate total number of kmer
    return df / s # calculate Kmer frequency by dividing count with total number of kmer 


# ## Functions for IntervalTree Package

# In[1]:


def INTERVALTREE_Visualize_Intervals( arr_intervals, sort_by_start = True, marker = 'o-' ) :
    ''' Visualize intervals (an object in IntervalTree python package) '''
    if type( arr_intervals ) is not np.ndarray : arr_intervals = np.array( arr_intervals, dtype = object )
    if sort_by_start : arr_intervals = arr_intervals[ np.argsort( list( interval[ 0 ] for interval in arr_intervals ) ) ]
    for index, interval in enumerate( arr_intervals ) :
        plt.plot( [ interval[ 0 ], interval[ 1 ] ], [ index, index ], marker, label = interval[ 2 ] )
    plt.legend( )


# ## Functions for Genomic Annotation Visualization

# In[1]:


def GTF_Visualize_Intervals_with_Exons_and_Introns( start, end, arr_intervals_of_interest = None, arr_pileup_count = None, figsize = ( 22, 6 ), show_legend = True, show_pileup_count_in_log_scale = False ) :
    ''' Visualize intervals (an object in IntervalTree python package) along with components (exons and introns) in 'it_comp' object and repeats (repeatmasker annotation) in 'it_re' object '''
    if arr_pileup_count is None : fig, ax_annotation = plt.subplots( 1, 1, figsize = figsize )
    else :
        fig, axes = plt.subplots( 2, 1, figsize = figsize, sharex = True )
        ax_pileup, ax_annotation = axes
        index, l_legend_line, l_legend_label = 0, list( ), list( )
        cycler_color = CYCLER( [ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf' ] )
    arr_repeats = np.array( list( it_re.overlap( start, end - 1 ) ), dtype = object ) 
    if len( arr_repeats ) > 0 : # if there is valid repeat element to visualize in a given window
        set_type_repeat = set( )
        for interval in arr_repeats :
            id_repeat, family_repeat, type_repeat, strand_repeat = interval[ 2 ]
            if type_repeat == 'SINE' :        alpha, lw, marker, color = 0.7, 4.5, 'x-', 'red' # graphic setting for each repeat element
            elif type_repeat == 'ERV' :       alpha, lw, marker, color = 0.7, 4.5, 'x-', 'orange'
            elif type_repeat == 'LINE' :      alpha, lw, marker, color = 0.7, 4.5, 'x-', 'green'
            elif type_repeat == 'DNA' :       alpha, lw, marker, color = 0.7, 4.5, 'x-', 'blue'
            elif type_repeat == 'Satellite' : alpha, lw, marker, color = 0.7, 4.5, 'x-', 'purple'
            else :                            alpha, lw, marker, color = 0.7, 4.5, 'x-', 'grey'
            if type_repeat not in set_type_repeat : # if current type_repeat does not have legend, add annotation to the legend
                l_legend_line.append( mpl.lines.Line2D( [ 0 ], [ 0 ], marker = 'x', color = color, lw = lw ) ), l_legend_label.append( type_repeat ) # add label to legend for annotation
                set_type_repeat.add( type_repeat )
            int_length_interval = interval[ 1 ] - interval[ 0 ]
            direction_arrow = 1 if strand_repeat == '+' else -1 
            ax_annotation.plot( [ interval[ 0 ], interval[ 1 ] ], [ index, index ], marker, lw = lw, color = color, alpha = alpha )
            ax_annotation.arrow( x = ( interval[ 0 ] + interval[ 1 ] ) / 2 - direction_arrow * int_length_interval / 10, y = index, dx = direction_arrow * int_length_interval / 5, dy = 0, fc = color, ec = color, length_includes_head = True, head_width = 0.8, head_length = int_length_interval / 5, shape = 'full', alpha = alpha )
        index += 1
    arr_intervals = np.array( list( it_comp.overlap( start, end - 1 ) ), dtype = object )
    if len( arr_intervals ) > 0 : # if there is valid component to visualize in a given window
        set_tx_id = set( )
        for interval in arr_intervals :
            for comp_info in dict_comp_id_to_l_gene_id__tx_id__type__strand[ interval[ 2 ] ] : set_tx_id.add( comp_info[ 1 ] )
        df_gene_tx = pd.DataFrame( list( [ dict_tx_id_to_gene_id[ tx_id ], tx_id ] for tx_id in set_tx_id ), columns = [ 'gene_id', 'tx_id' ] ).sort_values( [ 'gene_id', 'tx_id' ] )
        l_tx_id, l_gene_id = df_gene_tx.tx_id.values, df_gene_tx.gene_id.values
        l_l_interval_for_a_tx = list( list( [ dict_comp_id_to_start[ comp_id ], dict_comp_id_to_end[ comp_id ], comp_id, comp_type, dict_tx_id_to_strand[ tx_id ] ] for comp_id, comp_type in dict_tx_id_to_l_comp_id__type[ tx_id ] ) for tx_id in l_tx_id )
        gene_id_previous, color = l_gene_id[ 0 ], next( cycler_color )
        l_legend_line.append( mpl.lines.Line2D( [ 0 ], [ 0 ], color = color, lw = 3 ) ), l_legend_label.append( gene_id_previous )
        for l_interval_for_a_tx, tx_id, gene_id in zip( l_l_interval_for_a_tx, l_tx_id, l_gene_id ) :
            if gene_id != gene_id_previous : 
                color = next( cycler_color ) # change the color if gene_id changes
                l_legend_line.append( mpl.lines.Line2D( [ 0 ], [ 0 ], color = color, lw = 3 ) ), l_legend_label.append( gene_id ) # add label to legend for annotation
            gene_id_previous = gene_id
            arr_interval_for_a_tx = np.array( l_interval_for_a_tx, dtype = object )
            arr_interval_for_a_tx = arr_interval_for_a_tx[ arr_interval_for_a_tx[ :, 0 ].argsort( ) ] # sort intervals by its start position
            for interval in arr_interval_for_a_tx :
                if interval[ 3 ] == 'exon' : alpha, lw, marker = 1, 6, '-' # graphic setting for exon and introns
                else : alpha, lw, marker = 0.5, 3, '--'
                int_length_interval = interval[ 1 ] - interval[ 0 ]
                direction_arrow = 1 if interval[ 4 ] == '+' else -1 
                ax_annotation.plot( [ interval[ 0 ], interval[ 1 ] ], [ index, index ], marker, lw = lw, color = color, alpha = alpha )
                ax_annotation.arrow( x = ( interval[ 0 ] + interval[ 1 ] ) / 2 - direction_arrow * int_length_interval / 10, y = index, dx = direction_arrow * int_length_interval / 5, dy = 0, fc = color, ec = color, length_includes_head = True, head_width = 0.8, head_length = int_length_interval / 5, shape = 'full', alpha = alpha )
            index += 1
    if arr_intervals_of_interest is not None : # visualize given 'arr_intervals_of_interest'
        for interval in arr_intervals_of_interest :
            alpha, lw, marker, color = 1, 6, '-', next( cycler_color ) # graphic setting for intervals_of_interest
            l_legend_line.append( mpl.lines.Line2D( [ 0 ], [ 0 ], color = color, lw = 3 ) ), l_legend_label.append( interval[ 2 ] )
            int_length_interval = interval[ 1 ] - interval[ 0 ]
            direction_arrow = 1 if interval[ 3 ] == '+' else -1 
            ax_annotation.plot( [ interval[ 0 ], interval[ 1 ] ], [ index, index ], marker, lw = lw, color = color, alpha = alpha )
            ax_annotation.arrow( x = ( interval[ 0 ] + interval[ 1 ] ) / 2 - direction_arrow * int_length_interval / 10, y = index, dx = direction_arrow * int_length_interval / 5, dy = 0, fc = color, ec = color, length_includes_head = True, head_width = 0.8, head_length = int_length_interval / 5, shape = 'full', alpha = alpha )
    ax_annotation.set_ylim( -1, index + 1 )
    if show_legend : ax_annotation.legend( l_legend_line, l_legend_label )
            
    if arr_pileup_count is not None : 
        n_bases = end - start + 1
        int_sampling_a_base_every_n_bases = max( 1, int( n_bases / 10000 ) )
        x_range, data = np.arange( start, end + 1, int_sampling_a_base_every_n_bases ), arr_pileup_count[ start - 1 : end : int_sampling_a_base_every_n_bases ]
        x_axis = np.zeros( len( x_range ) )
        ax_pileup.plot( x_range, data, color = 'black' )
        ax_pileup.fill_between( x_range, data, x_axis, facecolor = 'black', interpolate = True, alpha = 0.5 )
        ax_pileup.set_xlim( start, end + 1 )
        if show_pileup_count_in_log_scale : ax_pileup.set_yscale( 'log' )


# ## Functions for Alignment Programs

# In[ ]:


def STAR_Read_Log( dir_glob ) :
    """read STAR final out logs file using a given glob string and return a dataframe"""
    dict_filename_to_dict_STAR_log = dict( )
    for dir_file in glob.glob(dir_glob):
        with open(dir_file) as file :
            l_line = file.read().split("\n")  
        dict_STAR_log = dict()
        for key,value in list(line.strip().split(" |\t") for line in Search_list_of_strings_with_multiple_query(l_line, "\t")):
            if '%' in value :
                value = value[:-1]
            try :
                value = float(value)
            except ValueError :
                pass
            dict_STAR_log[key] = value
        dict_filename_to_dict_STAR_log[dir_file.rsplit("/",1)[1]] = dict_STAR_log
    df_STAR_log = pd.DataFrame(dict_filename_to_dict_STAR_log)
    return df_STAR_log


# ## Functions for Annotations 

# In[ ]:


def ANNO_Ensembl_to_Entrez( df ) :
    ''' convert Ensembl gene_id to Entrez gene_id '''
    MAP.dict_a2b = df_anno_Ensembl_Entrez.Entrez_Gene_ID.to_dict( )
    df.index.name = 'gene_id'
    df = df.reset_index( )
    df[ 'gene_id' ] = df.gene_id.apply( MAP.Map_a2b )
    df = df.dropna( subset = [ 'gene_id' ] )
    df = df.set_index( 'gene_id' )
    return df


# ## Functions for STAR aligners

# In[ ]:


def STAR_Parse_final_out( dir_file, return_numeric = False ) :
    ''' parse final.out output file from STAR and return a series containing parsed values '''
    with open( dir_file ) as file : l_value = list( line.strip( ).split( ' |\t' ) for line in file.read( ).split( '\n' ) if ' |\t' in line )
    s = pd.DataFrame( l_value, columns = [ 'field', 'value' ] ).set_index( 'field' ).value
    if return_numeric :
        s.iloc[ : ] = list( value if ':' in value else float( value.replace( '%', '' ) ) for value in s.values )
    return s


# ## Functions for Fastq file

# In[ ]:


def FASTQ_Iterate( dir_file, return_only_at_index = 1 ) :
    """ # 2020-12-09 22:22:34 
    iterate through a given fastq file """
    if return_only_at_index is not None : return_only_at_index = return_only_at_index % 4 # 'return_only_at_index' value should be a value between 0 and 3
    bool_flag_file_gzipped = '.gz' in dir_file[ - 3 : ] # set a flag indicating whether a file has been gzipped.
    with gzip.open( dir_file, 'rb' ) if bool_flag_file_gzipped else open( dir_file ) as file :
        while True :
            record = [ file.readline( ).decode( )[ : -1 ] for index in range( 4 ) ] if bool_flag_file_gzipped else [ file.readline( )[ : -1 ] for index in range( 4 ) ]
            if len( record[ 0 ] ) == 0 : break
            if return_only_at_index is not None : yield record[ return_only_at_index ]
            else : yield record
                    
def FASTQ_Read( dir_file, return_only_at_index = None, return_generator = False ) : # 2020-08-18 22:31:31 
    ''' read a given fastq file into list of sequences or a dataframe (gzipped fastq file supported). 'return_only_at_index' is a value between 0 and 3 (0 = readname, 1 = seq, ...)
    'return_generator' : if True, return a generator that return a tuple of 4 length (representing a record with sequencing quality) or a value at the index given by "return_only_at_index".  '''
    if return_only_at_index is not None : return_only_at_index = return_only_at_index % 4 # 'return_only_at_index' value should be a value between 0 and 3
    bool_flag_file_gzipped = '.gz' in dir_file[ - 3 : ] # set a flag indicating whether a file has been gzipped.
    l_seq = list( )
    l_l_values = list( )
    file = gzip.open( dir_file, 'rb' ) if bool_flag_file_gzipped else open( dir_file )
    while True :
        record = [ file.readline( ).decode( )[ : -1 ] for index in range( 4 ) ] if bool_flag_file_gzipped else [ file.readline( )[ : -1 ] for index in range( 4 ) ]
        if len( record[ 0 ] ) == 0 : break
        if return_only_at_index is not None : l_seq.append( record[ return_only_at_index ] )
        else : l_l_values.append( [ record[ 0 ], record[ 1 ], record[ 3 ] ] )  
    file.close( )
    return l_seq if return_only_at_index is not None else pd.DataFrame( l_l_values, columns = [ 'readname', 'seq', 'quality' ] )


# In[ ]:


def FASTQ_Write( dir_file, dict_seq = None, dict_quality = None, df_fastq = None ) : 
    '''
    # 2021-03-07 15:53:19 
    Write FASTQ file with given dataframe or a pair of dict_seq and dict_quality (keys of both dictionaries should be the same) '''
    if dict_seq is None or dict_quality is None : 
        if df_fastq is None :
            print( 'required inputs are not given, exiting' )
            return -1
        else : # retrieve dictionary of sequences and quality from the dataframe
            df = df_fastq.set_index( 'readname' )
            dict_seq = df.seq.to_dict( )
            dict_quality = df.quality.to_dict( )
    flag_gzipped = dir_file.rsplit( '.', 1 )[ 1 ] == 'gz' # identify gzipped status
    with gzip.open( dir_file, 'wb' ) if flag_gzipped else open( dir_file, 'w' ) as file :
        for key in dict_seq :
            str_record = '\n'.join( [ '@' + key if key[ 0 ] != '@' else key, dict_seq[ key ], '+', dict_quality[ key ] ] ) + '\n'
            file.write( str_record.encode( ) if flag_gzipped else str_record )


# In[ ]:


def FASTQ_Split_with_adaptor( dir_file, dir_file_new, str_seq_adaptor ) :
    ''' # 2020-12-09 03:09:54 
    Split sequence record in the FASTQ file with a given adaptor sequence 
    '''
    flag_newfile_gzipped = dir_file_new.rsplit( '.', 1 )[ 1 ] == 'gz' # retrieve flag for gzipped status of the output file
    newfile = gzip.open( dir_file_new, 'wb' ) if flag_newfile_gzipped else open( dir_file_new, 'w' )
    for record in FASTQ_Read( dir_file, return_generator = True, return_only_at_index = None ) :
        name_read = record[ 0 ].split( ' ', 1 )[ 0 ]
        str_seq, str_qual = record[ 1 ], record[ 3 ]
        l_pos = [ 0 ] + STR.Find_all( record[ 1 ], str_seq_adaptor ) + [ len( record[ 1 ] ) ]
        for index in range( len( l_pos ) - 1 ) :
            str_record = '\n'.join( [ name_read + '_' + str( index ), str_seq[ l_pos[ index ] : l_pos[ index + 1 ] ], '+', str_qual[ l_pos[ index ] : l_pos[ index + 1 ] ] ] ) + '\n'
            newfile.write( str_record.encode( ) if flag_newfile_gzipped else str_record )
    newfile.close( )


# ## Functions for XML file

# In[ ]:


def XML_Read( dir_file, encoding = 'utf-8', character_for_replacing_invalid_characters = ' ' ) : # 2020-08-18 18:07:59 
    ''' read XML file with a given encoding as an ordered dictionary '''
    with open( dir_file, 'r', encoding = encoding ) as file :
        data = file.read( )
        try : dict_xml = xmltodict.parse( data )
        except ExpatError : 
            for int_char in np.arange( 32 ) :
                if int_char in [ 9, 10, 13 ] : continue # ignore white space characters
                data = data.replace( chr( int_char ), character_for_replacing_invalid_characters ) # remove special characters from xml content before parsing
            dict_xml = xmltodict.parse( data )
    return dict_xml


# ## Functions for Reading HMMER Output Files

# In[ ]:


def HMMER_Map_Positions_Between_Pairwise_Alignment( query_alignment, target_alignment, query_start, target_start, char_gap = '-' ) : # 2020-07-03 13:56:35 
    ''' return a dictionary that map positions (1-based coordinate) between the target sequence to the query sequence by using a pair-wise alignment result
    query_start, target_start: 1-based coordinates
    char_gap: a character representing a gap between a pair-wise alignment '''
    dict_target_pos_to_query_pos = dict( )
    int_start_query = 0 # position of a query sequence in the pair-wise alignment (1-based coordinate)
    int_start_target = 0 # position of a target sequence in the pair-wise alignment (1-based coordinate)
    for residue_query_alignment, residue_target_alignment in zip( query_alignment, target_alignment ) :
        if residue_query_alignment != char_gap and residue_target_alignment != char_gap : # no gap in either target or query sequences
            int_start_query += 1
            int_start_target += 1
            dict_target_pos_to_query_pos[ target_start + int_start_target - 1 ] = query_start + int_start_query - 1
        elif residue_target_alignment == char_gap : # gap in target sequence
            int_start_query += 1
        elif residue_query_alignment == char_gap : # gap in query sequence
            int_start_target += 1
            dict_target_pos_to_query_pos[ target_start + int_start_target - 1 ] = query_start + int_start_query - 1 # position in target sequence inside insertion will be mapped to the position where insertion begins
    return dict_target_pos_to_query_pos


# In[ ]:


def HMMER_HMMSEARCH_Read_output( dir_file ) :
    ''' read hmmsearch output file (output file of hmmsearch when jackhmmer-produced HMM files are used)  '''
    l_l_values = list( ) # list-based 2D array of values that will store alignment outputs for all hmmsearch output records
    with open( dir_file ) as file : 
        while True : # iterate line in the file
            line = file.readline( )
            if len( line ) == 0 : break
            else : content = line[ : -1 ]
            if len( content ) == 0 : continue # skip black line
            elif '#' == content[ 0 ] : continue # skip comments
            elif 'Query:' == content[ : 6 ] : acc_query = content.split( )[ 1 ] # retrive accession of the current query profile (human sequence's uniprot accession)
            elif '>>' == content[ : 2 ] : # for each alignment record
                acc_target = content[ 2: ].split( )[ 0 ] # retrive an identifier for the target protein (identifier is the first entry of white-space splitted fasta header of the target protein)
                line_header = file.readline( )
                if ' [No individual domains that' in line_header : continue  # if current alignment records does not contain domain-level alignment, continue on next record.
                file.readline( ) # discard header lines  # '   score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc', ' ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----',
                l_l_values_of_a_record = list( ) # list-based 2D array of values that will store an alignment output for a hmmsearch output record of a given query profile and a target protein
                while True :
                    line = file.readline( )
                    if len( line ) == 1 : break # if line contain only newline character, finish retriving alignment records in tabular format
                    else : content = line[ : -1 ]
                    l_l_values_of_a_record.append( [ acc_query, acc_target ] + content.split( ) )
                file.readline( ) # discard one line #  'Alignments for each domain:'
                for index in range( len( l_l_values_of_a_record ) ) :
                    cc = file.readline( ) # discard one line # == domain 1  score: 9.0 bits;  conditional E-value: 0.019', 
                    line_alignment = file.readline( ) # '                                     xxxxxxxxxxxxxxxxxxxxxxxxxxxxx....xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx RF' for jackhmmer HMM model after first iteration. HMM model from the first iteration does not contain this field, and it should be skipped.
                    str_query_alignment = line_alignment.split( )[ 2 ] if acc_query in line_alignment else file.readline( ).split( )[ 2 ] # '           sp|P10155|RO60_HUMAN-i4 372 vllvDvSGSMsasvsgsrltaakaaaalaaallrrgdrvglvaFdtr.vvdlpltprdsvlallaalaalpgGgTdigaalraalrrkar..rdvvvliTDgeendgpehptealae 485',
                    dd = file.readline( ) # discard one line # '                                       ++l+DvS S ++       + a+ +++ a+  l++g ++g++aF+ +  vd p++++  +  + a ++     g d+  alr+a +  +   +  v +i+Dg +  ++++  +++a+',
                    str_target_alignment = file.readline( ).split( )[ 2 ]
                    file.readline( ), file.readline( ) # discard two lines # '                                       799********33.....36899**********************88345656666544555543333...58*************888877778*******998888765.44443 PP', ''
                    l_l_values_of_a_record[ index ].extend( [ str_query_alignment, str_target_alignment ] )
                l_l_values.extend( l_l_values_of_a_record ) # add alignment information for current record to the 2D array
    df = pd.DataFrame( l_l_values ) # build a dataframe of alignment records.
    df.columns = [ 'query_accession', 'target_accession', 'index_match', 'significance', 'score', 'bias', 'conditional_Evalue', 'independent_Evalue', 'query_start', 'query_end', 'query_status', 'target_start', 'target_end', 'target_status', 'target_envelope_start', 'target_envelope_end', 'target_envelope_status', 'accuracy', 'query_alignment', 'target_alignment' ]
    l_col_interger_datatype = [ 'index_match', 'query_start', 'query_end', 'target_start', 'target_end', 'target_envelope_start', 'target_envelope_end' ] # column names with float datatype
    l_col_float_datatype = [ 'score', 'bias', 'conditional_Evalue', 'independent_Evalue', 'accuracy' ] # column names with integer datatype
    df.loc[ :, l_col_float_datatype ] = df.loc[ :, l_col_float_datatype ].astype( float ) # convert values in columns for numeric datatype to float datatype
    df.loc[ :, l_col_interger_datatype ] = df.loc[ :, l_col_interger_datatype ].astype( int ) # convert values in columns for numeric datatype to integer datatype
    return df


# In[ ]:


def HMMER_Add_CIGAR_String( df, dict_hmm_model_name_to_consensus_length ) : # 2020-05-20 17:02:41 
    ''' (Internal function) Convert HMMER alignment result ( 'seq_consensus', 'seq_target' columns ) to CIGAR string. Additionally, count number of matches, mismatches, insertion, and deletion in the alignment '''
    df.seq_consensus = df.seq_consensus.str.upper( )
    df.seq_target = df.seq_target.str.upper( )
    l_l_values = list( )
    for model_name, int_hmmfrom, int_hmm_to, str_consensus_alignment, str_target_alignment in df[ [ 'model_name', 'hmmfrom', 'hmm_to', 'seq_consensus', 'seq_target' ] ].values :
        len_seq_consensus = dict_hmm_model_name_to_consensus_length[ model_name ]
        int_length_seq = len( str_consensus_alignment )
        l_cigar_comp = list( ) # list of strings that will yield cigar string when joined 
        l_cigar_comp.append( [ 0, 'START' ] ) # put a component indicating the start of cigar string
        if int_hmmfrom > 1 : l_cigar_comp.append( [ ( int_hmmfrom - 1 ), "S" ] ) # introduce soft-cliping record to cigar string where concensus sequence was not mapped 
        bool_flag_insertion, bool_flag_deletion, bool_flag_match, bool_flag_mismatch = False, False, False, False
        for index in range( int_length_seq ) :
            char_consensus, char_target = str_consensus_alignment[ index ], str_target_alignment[ index ]
            str_operation = l_cigar_comp[ - 1 ][ 1 ]
            if char_consensus == '.' : # when deletion (respect to genome, deletion in the consensus sequence) occured at current position
                if str_operation == 'D' : l_cigar_comp[ - 1 ][ 0 ] += 1
                else : l_cigar_comp.append( [ 1, 'D' ] )
            elif char_target == '-' : # when insertion (respect to genome, insertion in the consensus sequence) occured at current position
                if str_operation == 'I' : l_cigar_comp[ - 1 ][ 0 ] += 1
                else : l_cigar_comp.append( [ 1, 'I' ] )
            elif char_consensus == char_target : # when match occured at current position
                if str_operation == '=' : l_cigar_comp[ - 1 ][ 0 ] += 1
                else : l_cigar_comp.append( [ 1, '=' ] )
            else : # when mismatch occured at current position
                if str_operation == 'X' : l_cigar_comp[ - 1 ][ 0 ] += 1
                else : l_cigar_comp.append( [ 1, 'X' ] )
        #     else : # when match or mismatch occured at current position
        #         if str_operation == 'M' : l_cigar_comp[ - 1 ][ 0 ] += 1
        #         else : l_cigar_comp.append( [ 1, 'M' ] )
        if int_hmm_to < len_seq_consensus : l_cigar_comp.append( [ len_seq_consensus - int_hmm_to, "S" ] ) # # introduce soft-cliping record to cigar string where concensus sequence was not mapped 
        l_cigar_comp = l_cigar_comp[ 1 : ] # exclude the component that was used to indicate the start of cigar string
        arr_cigar_comp = np.array( l_cigar_comp, dtype = object )
        int_count_matched = arr_cigar_comp[ :, 0 ][ arr_cigar_comp[ :, 1 ] == '=' ].sum( )
        int_count_mismatched = arr_cigar_comp[ :, 0 ][ arr_cigar_comp[ :, 1 ] == 'X' ].sum( )
        int_count_insertion = arr_cigar_comp[ :, 0 ][ arr_cigar_comp[ :, 1 ] == 'I' ].sum( )
        int_count_deletion = arr_cigar_comp[ :, 0 ][ arr_cigar_comp[ :, 1 ] == 'D' ].sum( )
        str_cigar = ''.join( list( str( entry[ 0 ] ) + entry[ 1 ] for entry in l_cigar_comp ) )
        l_l_values.append( [ str_cigar, int_count_matched, int_count_mismatched, int_count_insertion, int_count_deletion ] )
    df = df.join( pd.DataFrame( l_l_values, columns = [ 'CIGAR_String', 'n_match', 'n_mismatch', 'n_insertion', 'n_deletion' ] ) )
    return df


# In[ ]:


def HMMER_Read_nhmmscan_output_with_alignment( dir_file ) : # 2020-05-20 17:02:41 
    ''' Read nhmmscan output file with alignment information. read only search summary, not alignments, and return a dataframe containing the summary '''
    with open( dir_file ) as file :
        l_l_value__search_result = list( )
        l_l_value__alignment_result = list( )
        while True : # retrive a directory to a file containing HMM database in the header
            line = file.readline( )
            if len( line ) == 0 : break 
            else : content = line[ : -1 ] 
            if '# target HMM database:' in content : # check whether current line containing HMM db directory
                dir_file_hmm_db = content.split( '# target HMM database: ' )[ 1 ].strip( )
                dict_fasta = FASTA_Read( os.popen( 'hmmemit -C {}'.format( dir_file_hmm_db ) ) )
                dict_hmm_model_name_to_consensus_length = dict( ( hmm_model_name.rsplit( '-consensus', 1 )[ 0 ], len( dict_fasta[ hmm_model_name ] ) ) for hmm_model_name in dict_fasta ) # retrive length of consensus sequences in a given HMM db, while removing '-consensus' tag from the hmm_model_name
                break
        while True :
            while True :
                line = file.readline( )
                if len( line ) == 0 : break 
                else : content = line[ : -1 ] 
                if len( content ) > 0 and content[ 0 ] != "#" : # ignore empty content and comments
                    if 'Query:' in content : 
                        str_query_sequence_name = content.split( 'Query:' )[ 1 ].split( '[L=' )[ 0 ].strip( )
                        int_query_sequence_length = int( content.split( '[L=' )[ 1 ].split( ']', 1 )[ 0 ] )
                        bool_flag_new_query_sequence_start = True
                    elif '-------' in content and bool_flag_new_query_sequence_start : 
                        str_field_locations = content
                        l_l_field_annotation_start_end = STR.Read_Field_Locations( str_field_locations )
                        slice_field_e_value = slice( l_l_field_annotation_start_end[ 0 ][ 0 ], l_l_field_annotation_start_end[ 0 ][ 1 ] )
                        slice_field_score = slice( l_l_field_annotation_start_end[ 1 ][ 0 ], l_l_field_annotation_start_end[ 1 ][ 1 ] )
                        slice_field_bias = slice( l_l_field_annotation_start_end[ 2 ][ 0 ], l_l_field_annotation_start_end[ 2 ][ 1 ] )
                        int_start_position__field_start = l_l_field_annotation_start_end[ 4 ][ 1 ] - len( str( int_query_sequence_length ) ) 
                        slice_field_start = slice( int_start_position__field_start, l_l_field_annotation_start_end[ 4 ][ 1 ] )
                        slice_field_end = slice( l_l_field_annotation_start_end[ 5 ][ 1 ] - len( str( int_query_sequence_length ) ), l_l_field_annotation_start_end[ 5 ][ 1 ] )
                        slice_field_model = slice( l_l_field_annotation_start_end[ 3 ][ 0 ], int_start_position__field_start - 1 )
                        bool_flag_new_query_sequence_start = False
                        break
                    elif '>>' in content : # retrive alignment records
                        bool_flag_new_alignment_start = True 
                        str_aligned_model_name = content[ 3 : ].strip( )
                    elif '-------' in content and bool_flag_new_alignment_start :
                        str_aligned_section_field_locations = content
                        l_l_aligned_section_field_annotation_start_end = STR.Read_Field_Locations( str_aligned_section_field_locations )
                        content = file.readline( )[ : -1 ]
                        bool_flag_new_alignment_start = False
                        file.readline( )
                        if 'Alignment:' not in file.readline( ) : print( 'error' )
                        if 'score:' not in file.readline( ) : print( 'error' )
                        content_consensus = file.readline( )[ : -1 ]
                        _                 = file.readline( )[ : -1 ]
                        content_target    = file.readline( )[ : -1 ]
                        content_quality   = file.readline( )[ : -1 ]
                        index_start, value_counter = 0, 0 # retrive sequence start index and end index
                        for value in content_consensus.split( ' ' )[ 1 : ] :
                            if value_counter < 2 :
                                if len( value ) > 0 : value_counter += 1
                                index_start += ( 1 + len( value ) )
                            else :
                                index_start += 1
                                index_end = index_start + len( value ) 
                                break
                        l_l_value__alignment_result.append( [ str_query_sequence_name, str_aligned_model_name, content[ l_l_aligned_section_field_annotation_start_end[ 0 ][ 0 ] - 2 ], float( content[ l_l_aligned_section_field_annotation_start_end[ 0 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 0 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 1 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 1 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 2 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 2 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 3 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 3 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 4 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 4 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 4 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 5 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 5 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 5 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 6 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 6 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 6 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 7 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 7 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 7 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 8 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 8 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 8 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 9 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 9 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 9 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 10 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 10 ][ 1 ] ] ), content_consensus[ index_start : index_end ], content_target[ index_start : index_end ], content_quality[ index_start : index_end ] ] )
            while True :
                line = file.readline( )
                if len( line ) == 0 : break 
                else : content = line[ : -1 ]
                if len( content ) == 0 : break # finished the loop when empty content is detected
                if '------ inclusion threshold ------' in content : continue # ignore the content indicating inclusion threshold
                l_l_value__search_result.append( [ float( content[ slice_field_e_value ] ), float( content[ slice_field_score ] ), float( content[ slice_field_bias ] ), int( content[ slice_field_start ] ), int( content[ slice_field_end ] ), str_query_sequence_name, content[ slice_field_model ].strip( ) ] )
            if len( line ) == 0 : break 
    df_search_result = pd.DataFrame( l_l_value__search_result, columns = [ 'e_value', 'score', 'bias', 'start', 'end', 'query_sequence_name', 'model_name' ] )
    df_alignment_result = HMMER_Add_CIGAR_String( pd.DataFrame( l_l_value__alignment_result, columns = [ 'query_sequence_name', 'model_name', 'inclusion', 'score', 'bias', 'Evalue', 'hmmfrom', 'hmm_to', 'hmm_info', 'alifrom', 'ali_to', 'ali_info', 'envfrom', 'env_to', 'env_info', 'mod_len', 'acc', 'seq_consensus', 'seq_target', 'seq_quality' ] ), dict_hmm_model_name_to_consensus_length )
    return df_search_result, df_alignment_result


# In[ ]:


def HMMER_Read_nhmmer_output_with_alignment( dir_file ) : 
    ''' # 2021-02-04 11:20:08 
    Read nhmmer output file with alignment information. read search summary, alignments, and return two dataframes containing the summary and alignment result '''
    bool_flag_new_alignment_start = False # in case that there is no alignment
    with open( dir_file ) as file :
        l_l_value__search_result = list( )
        l_l_value__alignment_result = list( )
        while True : # retrive a directory to a file containing HMM database in the header
            line = file.readline( )
            if len( line ) == 0 : break 
            else : content = line[ : -1 ] 
            if '# query file:' in content : # check whether current line containing HMM db directory
                dir_file_hmm_db = content.split( '# query file: ' )[ 1 ].strip( )
                if dir_file_hmm_db.rsplit( '.', 1 )[ 1 ].lower( ) in [ 'fa', 'fasta' ] : # if query file is a fasta file
                    dict_fasta = FASTA_Read( dir_file_hmm_db )
                    dict_hmm_model_name_to_consensus_length = dict( ( dir_file_hmm_db.rsplit( '/', 1 )[ 1 ].rsplit( '.', 1 )[ 0 ], len( dict_fasta[ hmm_model_name ] ) ) for hmm_model_name in dict_fasta ) # retrive length of sequence # model name is filename
                else : # if query file is hmm file
                    dict_fasta = FASTA_Read( os.popen( 'hmmemit -C {}'.format( dir_file_hmm_db ) ) )
                    dict_hmm_model_name_to_consensus_length = dict( ( hmm_model_name.rsplit( '-consensus', 1 )[ 0 ], len( dict_fasta[ hmm_model_name ] ) ) for hmm_model_name in dict_fasta ) # retrive length of consensus sequences in a given HMM db, while removing '-consensus' tag from the hmm_model_name
                break
        while True :
            while True :
                line = file.readline( )
                if len( line ) == 0 : break 
                else : content = line[ : -1 ] 
                if len( content ) > 0 and content[ 0 ] != "#" : # ignore empty content and comments
                    if 'Query:' in content : 
                        str_query_sequence_name = content.split( 'Query:' )[ 1 ].split( '[M=' )[ 0 ].strip( )
                        int_query_sequence_length = int( content.split( '[M=' )[ 1 ].split( ']', 1 )[ 0 ] )
                        bool_flag_search_result_table_start = True
                    elif '-------' in content and bool_flag_search_result_table_start : 
                        str_field_locations = content
                        l_l_field_annotation_start_end = STR.Read_Field_Locations( str_field_locations )
                        slice_field_e_value = slice( l_l_field_annotation_start_end[ 0 ][ 0 ], l_l_field_annotation_start_end[ 0 ][ 1 ] )
                        slice_field_score = slice( l_l_field_annotation_start_end[ 1 ][ 0 ], l_l_field_annotation_start_end[ 1 ][ 1 ] )
                        slice_field_bias = slice( l_l_field_annotation_start_end[ 2 ][ 0 ], l_l_field_annotation_start_end[ 2 ][ 1 ] )
                        slice_field_Sequence = slice( l_l_field_annotation_start_end[ 3 ][ 0 ], l_l_field_annotation_start_end[ 3 ][ 1 ] )
                        slice_field_start = slice( l_l_field_annotation_start_end[ 3 ][ 1 ] + 1, l_l_field_annotation_start_end[ 4 ][ 1 ] )
                        slice_field_end = slice( l_l_field_annotation_start_end[ 4 ][ 1 ] + 1, l_l_field_annotation_start_end[ 5 ][ 1 ] )
                        slice_field_Description = slice( l_l_field_annotation_start_end[ 5 ][ 1 ] + 1, None )
                        bool_flag_search_result_table_start = False
                        int_n_extra_characters_for_Sequence = 0 # number of extra characters for Sequence annotation
                        break
                    elif '>>' in content : # retrive alignment records
                        bool_flag_new_alignment_start = True 
                        int_description_start = content.find( 'dna:' )
                        str_Sequence, str_Description = content[ 3 : int_description_start ].strip( ), content[ int_description_start : ].strip( )
                    elif '-------' in content and bool_flag_new_alignment_start :
                        str_aligned_section_field_locations = content
                        l_l_aligned_section_field_annotation_start_end = STR.Read_Field_Locations( str_aligned_section_field_locations )
                        content = file.readline( )[ : -1 ]
                        bool_flag_new_alignment_start = False
                        file.readline( ) # skip blank and other lines that does not contain meaningful information
                        if 'Alignment:' not in file.readline( ) : print( 'error' )
                        if 'score:' not in file.readline( ) : print( 'error' )
                        str_first_line_alignment = file.readline( )[ : -1 ] # retrieve first line of the alignment
                        content_consensus = file.readline( )[ : -1 ] if ' RF' in str_first_line_alignment else str_first_line_alignment # if a mask of reference sequence showing indels is present in the first line of the alignment record, use the next line for consensus sequence
                        _                 = file.readline( )[ : -1 ]
                        content_target    = file.readline( )[ : -1 ]
                        content_quality   = file.readline( )[ : -1 ]
                        index_start, value_counter = 0, 0 # retrive sequence start index and end index
                        for value in content_consensus.split( ' ' )[ 1 : ] :
                            if value_counter < 2 :
                                if len( value ) > 0 : value_counter += 1
                                index_start += ( 1 + len( value ) )
                            else :
                                index_start += 1
                                index_end = index_start + len( value ) 
                                break
                        l_l_value__alignment_result.append( [ str_query_sequence_name, str_Sequence, str_Description, content[ l_l_aligned_section_field_annotation_start_end[ 0 ][ 0 ] - 2 ], float( content[ l_l_aligned_section_field_annotation_start_end[ 0 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 0 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 1 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 1 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 2 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 2 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 3 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 3 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 4 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 4 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 4 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 5 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 5 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 5 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 6 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 6 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 6 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 7 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 7 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 7 ][ 1 ] ] ), int( content[ l_l_aligned_section_field_annotation_start_end[ 8 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 8 ][ 1 ] ] ), content[ l_l_aligned_section_field_annotation_start_end[ 8 ][ 1 ] + 1 : l_l_aligned_section_field_annotation_start_end[ 9 ][ 0 ] - 1 ], int( content[ l_l_aligned_section_field_annotation_start_end[ 9 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 9 ][ 1 ] ] ), float( content[ l_l_aligned_section_field_annotation_start_end[ 10 ][ 0 ] : l_l_aligned_section_field_annotation_start_end[ 10 ][ 1 ] ] ), content_consensus[ index_start : index_end ], content_target[ index_start : index_end ], content_quality[ index_start : index_end ] ] )
            while True :
                line = file.readline( )
                if len( line ) == 0 : break 
                else : content = line[ : -1 ]
                if len( content ) == 0 : break # finished the loop when empty content is detected
                if '------ inclusion threshold ------' in content : continue # ignore the content indicating the inclusion threshold
                try : int( content[ slice_field_start ] )
                except :
                    while True :
                        int_n_extra_characters_for_Sequence += 1 # if more characters should be used to retrive 'Sequence' field, use one more character and check whether an entire 'Sequence' can be retrived.
                        if content[ l_l_field_annotation_start_end[ 3 ][ 1 ] + int_n_extra_characters_for_Sequence ] == ' ' : break
                    slice_field_Sequence = slice( l_l_field_annotation_start_end[ 3 ][ 0 ], l_l_field_annotation_start_end[ 3 ][ 1 ] + int_n_extra_characters_for_Sequence )
                    slice_field_start = slice( l_l_field_annotation_start_end[ 3 ][ 1 ] + 1 + int_n_extra_characters_for_Sequence, l_l_field_annotation_start_end[ 4 ][ 1 ] )
                l_l_value__search_result.append( [ float( content[ slice_field_e_value ] ), float( content[ slice_field_score ] ), float( content[ slice_field_bias ] ), content[ slice_field_Sequence ].strip( ), int( content[ slice_field_start ] ), int( content[ slice_field_end ] ), content[ slice_field_Description ].strip( ), str_query_sequence_name ] )

            if len( line ) == 0 : break 
    df_search_result = pd.DataFrame( l_l_value__search_result, columns = [ 'e_value', 'score', 'bias', 'Sequence', 'start', 'end', 'Description', 'model_name' ] )
    df_alignment_result = HMMER_Add_CIGAR_String( pd.DataFrame( l_l_value__alignment_result, columns = [ 'model_name', 'Sequence', 'Description', 'inclusion', 'score', 'bias', 'Evalue', 'hmmfrom', 'hmm_to', 'hmm_info', 'alifrom', 'ali_to', 'ali_info', 'envfrom', 'env_to', 'env_info', 'mod_len', 'acc', 'seq_consensus', 'seq_target', 'seq_quality' ] ), dict_hmm_model_name_to_consensus_length )
    return df_search_result, df_alignment_result


# In[ ]:


def NHMMER_Plot_Frequency( df_alignment_result, len_hmm, n_seq = None, thres_frequency_for_consensus_seq = 0.6 ) :
    """
    # 2021-02-04 14:16:44 
    plot frequency of nucleotide residue for each position of a HMM model used in the nhmmer
    'df_alignment_result' : one of the output dataframes of the 'HMMER_Read_nhmmer_output_with_alignment' function
    'len_hmm' : length of the HMM model
    'thres_frequency_for_consensus_seq' : threshold (float, 0 ~ 1) for base calling of consensus sequences from residue frequencies of the HMM model's positions
    """
    df = df_alignment_result
    if n_seq is None :
        n_seq = len( df ) # plot all aligned sequences by default
    str_residues = "ATGC-" # define list of valid nucleotide residues
    dict_residue_to_index = dict( ( residue, index ) for index, residue in enumerate( str_residues ) )
    arr_count = np.zeros( ( len( str_residues ), len_hmm ) )
    for record in df.iloc[ : n_seq ].to_dict( orient = 'record' ) :
        seq_hmm = record[ "seq_consensus" ]
        seq_target = record[ "seq_target" ]
        int_start_hmm = record[ "hmmfrom" ] - 1 - 1 # 0-based coordinate, before (-1) the first index
        int_start_target = - 1
        for residue_hmm, residue_target in zip( seq_hmm, seq_target ) :
            if residue_hmm == '.' :
                int_start_target += 1

            elif residue_target == '-' :
                int_start_hmm += 1
                arr_count[ dict_residue_to_index[ residue_target ], int_start_hmm ] += 1

            else :
                int_start_hmm += 1
                int_start_target += 1
                arr_count[ dict_residue_to_index[ residue_target ], int_start_hmm ] += 1

    arr_frequency = arr_count / n_seq # calculate residue frequency
    # plot residue frequencies
    fig, ax = plt.subplots( 1, 1, figsize = ( 15, 5 ) )
    for index, arr in enumerate( arr_frequency ) :
        ax.plot( arr, label = str_residues[ index ] )
    ax.legend( ) 
    
    arr_residue = np.full( len_hmm, 'N', dtype = object ) # 1d array containing consensus-base-called nucleotide residues
    for index_residue, index_pos in zip( * np.where( arr_frequency > thres_frequency_for_consensus_seq ) ) :
        arr_residue[ index_pos ] = str_residues[ index_residue ]
    str_consensus_seq = ''.join( arr_residue ) # retrieve consensus sequence
    return str_consensus_seq


# ## Functions for BLAST

# In[ ]:


def BLAST_Parse_BTOP_String( str_btop, query_seq = None, subject_seq = None ) : # 2020-07-02 17:48:53 
    ''' By using a Blast traceback operations (BTOP) String produced by BLASTp and query_seq or subject_seq, build aligned query and subject sequences. '''
    if query_seq is None and subject_seq is None : return -1 # either 'query_seq' or 'subject_seq' should be given to build aligned query and subject sequences.
    bool_flag_use_query_seq = query_seq is not None
    query_seq_aligned = ''
    subject_seq_aligned = ''
    str_number = ''
    index = 0
    iter_str_btop = iter( str_btop ) # make an iterator for the BTOP string.
    while True :
        try : char = next( iter_str_btop )
        except : break
        if char in '0123456789' : str_number += char # if encountered is a number character of a number indicating the number of local matches, add up these number characters
        else : 
            if len( str_number ) > 0 : # if there was a local match, retrieve a locally matched sequence from either query or subject sequences and append the sequence to aligned aligned query and subject sequences.
                n_matches = int( str_number )
                if bool_flag_use_query_seq : seq_matched = query_seq[ index : index + n_matches ]
                else : seq_matched = subject_seq[ index : index + n_matches ]
                index += n_matches # update index indicating position on eigher query or subject sequence.
                query_seq_aligned += seq_matched
                subject_seq_aligned += seq_matched
                str_number = ''
            char_2 = next( iter_str_btop ) # if mismatch or gap is encountered, update aligned sequences accordingly
            query_seq_aligned += char
            subject_seq_aligned += char_2
            if char == '-' and not bool_flag_use_query_seq : index += 1
            elif char_2 == '-' and bool_flag_use_query_seq : index += 1
            elif char != '-' and char_2 != '-' : index += 1
    if len( str_number ) > 0 : # if there was a local match just before the end of a given btop string, retrieve a locally matched sequence from either query or subject sequences and append the sequence to aligned aligned query and subject sequences.
        n_matches = int( str_number )
        if bool_flag_use_query_seq : seq_matched = query_seq[ index : index + n_matches ]
        else : seq_matched = subject_seq[ index : index + n_matches ]
        query_seq_aligned += seq_matched
        subject_seq_aligned += seq_matched
    return query_seq_aligned, subject_seq_aligned


# In[ ]:


def BLAST_Read( dir_file, dict_qaccver_to_seq = None, dict_saccver_to_seq = None, ** dict_Select ) : # 2020-11-15 17:21:00 
    """ Read BLASTp tabular output with following columns [ 'qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'btop' ]
    Subset BLASTp output with 'dict_Select' using PD_Select function before BTOP string parsing.
    For BTOP String parsing, either 'dict_qaccver_to_seq' or 'dict_saccver_to_seq' should be given as a template sequence. """
    if dict_qaccver_to_seq is None and dict_saccver_to_seq is None : return -1 # For BTOP String parsing, either 'dict_qaccver_to_seq' or 'dict_saccver_to_seq' should be given as a template sequence. 
    bool_flag_use_query_seq = dict_qaccver_to_seq is not None # set a flag for using query_seq
    try : 
        df_blastp = pd.read_csv( dir_file, sep = '\t', header = None )
    except pd.errors.EmptyDataError :
        print( "{dir_file} file contains empty blast output".format( dir_file = dir_file ) ); return -1 
    df_blastp.columns = [ 'qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'btop' ]
    df_blastp = PD_Select( df_blastp, ** dict_Select ) # retrieve search result of multiple-sequence-aligned sequences
    l_l_value = list( ) # build aligned query and subject sequences for blastp search results
    for qaccver, qstart, qend, saccver, sstart, send, btop in df_blastp[ [ 'qaccver', 'qstart', 'qend', 'saccver', 'sstart', 'send', 'btop' ] ].values :
        query_seq_aligned, subject_seq_aligned = BLAST_Parse_BTOP_String( btop, query_seq = dict_qaccver_to_seq[ qaccver ][ qstart - 1 : qend ] ) if bool_flag_use_query_seq else BLAST_Parse_BTOP_String( btop, subject_seq = dict_saccver_to_seq[ saccver ][ sstart - 1 : send ] ) # parse BTOP string according to the given type of sequence
        l_l_value.append( [ query_seq_aligned, subject_seq_aligned ] )
    df_blastp = df_blastp.reset_index( drop = True ).join( pd.DataFrame( l_l_value, columns = [ 'query_seq_aligned', 'subject_seq_aligned' ] ) )
    return df_blastp


# ## Functions for CIF format structure files

# In[ ]:


def CIF_Read( dir_file ) : # 2020-07-20 20:31:16 
    """ Read CIF format file (gzipped file supported) and return it as a dictionary containing dictionaries representing data_blocks, which contains category names as keys and DataFrames and Series as values. """
    bool_flag_gzipped = dir_file[ -3 : ] == '.gz' # flag for gzipped status of a given cif file
    dict_cif = dict( )
    str_spanning_multiple_line = None
    dict_data_block = dict( )
    l_l_value = list( )
    with gzip.open( dir_file, 'rb' ) if bool_flag_gzipped else open( dir_file ) as file :
        while True :
            line = file.readline( ).decode( ) if bool_flag_gzipped else file.readline( ) 
            if len( line ) == 0 : break
            else : line = line.strip( ) # read one line and strip empty characters at both ends
            if len( line ) == 0 : continue # ignore empty line
            if str_spanning_multiple_line is not None : # when the current line belongs to a multiple-lines-spanning string
                if line[ 0 ] == ';' : # when the end of a multiple-lines-spanning string is detected
                    if bool_flag_datetype_category_is_tabular and len( l_l_value[ 0 ] ) == len( l_l_value[ -1 ] ) : l_l_value.append( [ str_spanning_multiple_line ] ) # add a value spanning multiple lines as a first value of a new row if category datatype is tabular datatype and the number of values in the previous row is the same as the number of attribute names
                    else : l_l_value[ -1 ].append( str_spanning_multiple_line ) # add a value spanning multiple lines to the current row (eigher key-value or tabular category datatype)
                    str_spanning_multiple_line = None # flush a string containing value spanning multiple-lines
                else : str_spanning_multiple_line += line.replace( '\t', ' ' )  # when '#' character inside a value spanning multiple lines
            elif 'data_' == line[ : 5 ] : # when new data_block starts, add 'dict_data_block' to 'dict_cif' with 'str_name_data_block' and empties 'dict_data_block'
                if len( dict_data_block ) > 0 : dict_cif[ str_name_data_block ] = dict_data_block
                str_name_data_block = line[ 5 : ]
                dict_data_block = dict( )
            elif '#' == line[ 0 ] : # when new category start, add category data (Series format for key-value data and DataFrame format for tabular data) to 'dict_data_block' with 'str_name_category'
                if len( l_l_value ) > 0 :
                    str_category = '\n'.join( list( '\t'.join( l_value ) for l_value in l_l_value ) ) + '\n' # converting data to a 'file' and parse it using pandas.read_csv allows automatic data type detection and conversion when building pandas data
                    dict_data_block[ str_name_category ] = pd.read_csv( StringIO( str_category ), sep = '\t', low_memory = False, quoting = csv.QUOTE_NONE, keep_default_na = False ) if bool_flag_datetype_category_is_tabular else pd.read_csv( StringIO( str_category ), sep = '\t', squeeze = True, header = None, index_col = 0, low_memory = False, quoting = csv.QUOTE_NONE, keep_default_na = False ) # csv.QUOTE_NONE to deal with single quotation mark inside quotation # keep_default_na = False (prevent automatic conversion of NA chain to np.nan)
                bool_flag_datetype_category_is_tabular = False # 'bool_flag_datetype_category_is_tabular' is False by default (default category datatype is key-value datatype)
                l_l_value = list( )
            elif 'loop_' == line[ : 5 ] : # when 'loop_' directive is used, set 'bool_flag_datetype_category_is_tabular' to True (category datatype is tabular datatype)
                bool_flag_datetype_category_is_tabular = True
                l_l_value.append( [ ] ) # add a 'row' for collecting attribute names
            elif '_' == line[ 0 ] : # when CIF token name detected
                line = line[ 1 : ]
                if bool_flag_datetype_category_is_tabular : # if category datatype is tabular, retrieve attribute name
                    str_name_category, str_name_attribute = line.split( '.' )
                    l_l_value[ 0 ].append( str_name_attribute )
                else : # if category datatype is key-value, retrive attribute name and value (each 'row' in 'l_l_value' contains attribute name and value)
                    l_l_value.append( list( value.replace( '\t', ' ' ) for value in PARSE_Empty_Space_Delimited_with_Quotations( line ) ) ) # split a line into values using empty-space delimiter while ignoring delimiters inside quotations. # remove '\t' from value to be compativle with parsing method
                    str_name_category, l_l_value[ -1 ][ 0 ] = l_l_value[ -1 ][ 0 ].split( '.' ) # retrieve category name and attribute name from the current CIF token name
            elif line[ 0 ] == ';' : str_spanning_multiple_line = line[ 1 : ].replace( '\t', ' ' ) # when the start of a multiple-lines-spanning string is detected
            elif bool_flag_datetype_category_is_tabular : # for tabular category data
                l_value = list( value.replace( '\t', ' ' ) for value in PARSE_Empty_Space_Delimited_with_Quotations( line ) ) # split a line into values using empty-space delimiter while ignoring delimiters inside quotations. # remove '\t' from value to be compativle with parsing method
                if len( l_l_value[ 0 ] ) == len( l_l_value[ -1 ] ) : l_l_value.append( l_value ) # if the number of values in the previous row is the same as the number of attribute names, add the parsed row as a new row
                else : l_l_value[ -1 ].extend( l_value ) # if the number of values in the previous row is not equal to the number of attribute names, append the parsed row to the current row (these values of the same row are separated by a multiple-line-spanning value)
            else : l_l_value[ -1 ].append( line.replace( '\t', ' ' ) ) # for key-value category data, add a value spanning an entire line to the current row # remove '\t' from value to be compativle with parsing method
        if len( l_l_value ) > 0 : # when file ends, add category data to 'dict_data_block', and add data_block to 'dict_cif'
            str_category = '\n'.join( list( '\t'.join( l_value ) for l_value in l_l_value ) ) + '\n' # converting data to a 'file' and parse it using pandas.read_csv allows automatic data type detection and conversion when building pandas data
            dict_data_block[ str_name_category ] = pd.read_csv( StringIO( str_category ), sep = '\t', low_memory = False ) if bool_flag_datetype_category_is_tabular else pd.read_csv( StringIO( str_category ), sep = '\t', squeeze = True, header = None, index_col = 0, low_memory = False )
        if len( dict_data_block ) > 0 : dict_cif[ str_name_data_block ] = dict_data_block
    return dict_cif


# In[ ]:


def CIF_Write( dir_file, dict_cif ) : # 2020-07-13 23:05:05 
    """ write CIF file in the given 'dict_cif'. 
    *the output file closely resembles standard CIF file, but it does not have line length limit (no value spans multiple lines enclosed by semicolons ';'). This will prevent some CIF readers, such as mkdssp or chimera from reading the metadata (not 3D coordinate data). Therefore, saving only 'atom_site' category is recommended """
    with open( dir_file, 'w' ) as file :
        for str_name_data_block in dict_cif : # for each data block
            dict_data_block = dict_cif[ str_name_data_block ]
            file.write( 'data_' + str_name_data_block + '\n' )
            for str_name_category in dict_data_block : # for each category
                file.write( '#\n' )
                data = dict_data_block[ str_name_category ] # retrive category data (either Series or DataFrame)
                if isinstance( data, pd.Series ) :
                    arr = data.reset_index( ).values.astype( str ).astype( object )
                    arr[ :, 0 ] = '_' + str_name_category + '.' + arr[ :, 0 ]
                else :
                    file.write( 'loop_\n' )
                    file.write( '\n'.join( '_' + str_name_category + '.' + data.columns.values.astype( object ) ) + '\n' ) # write attributes first for tabular data
                    arr = data.round( 4 ).values.astype( str ).astype( object )
                str_format = ''.join( list( '{:<' + str( int_column_width ) + 's}' for int_column_width in np.array( list( list( len( value ) for value in arr_value ) for arr_value in arr ) ).max( axis = 0 ) + ( 3 if isinstance( data, pd.Series ) else 1 ) ) ) # retrive str_format by using maximum column_width for each column. minimum empty space for Series data is '   ', while that for tabular data is ' '
                file.write( '\n'.join( list( str_format.format( * arr_value ) for arr_value in arr ) ) + '\n' )
            file.write( '#\n' ) # write '#' at the end of the data_block


# In[ ]:


def CIF_from_PDB( df, l_label_entity_id = None ) : # 2020-08-08 15:21:53 
    """ Convert PDB dataframe (returned by PDB.Read_Single_Module) to mmCIF atom_site category dataframe. # default chain_id is 'A', if chain_id is empty in the PDB record """
    l_col_required = [ 'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'auth_asym_id' ] # necessary columns in a mmCIF format file
    df.rename( columns = { 'ATOM_or_HETATM' : 'group_PDB', 'Atom_serial_number' : 'id', 'Atom_name' : 'label_atom_id', 'Alternate_location_indicator' : 'label_alt_id', 'Residue_name' : 'label_comp_id', 'Chain_identifier' : 'label_asym_id', 'Residue_sequence_number' : 'label_seq_id', 'X_orthogonal_Å_coordinate' : 'Cartn_x', 'Y_orthogonal_Å_coordinate' : 'Cartn_y', 'Z_orthogonal_Å_coordinate' : 'Cartn_z' }, inplace = True )
    df.drop( columns = [ 'Code_for_insertions_of_residues' ], inplace = True )
    df[ 'type_symbol' ] = list( value[ 0 ] for value in df[ 'label_atom_id' ].values )
    df[ 'label_entity_id' ] = 1 if l_label_entity_id is None else l_label_entity_id # default value of l_label_entity_id is 1
    df.label_alt_id = df.label_alt_id.replace( '', '.' )
    df.label_asym_id = df.label_asym_id.replace( '', 'A' ) # default chain_id (label_asym_id) is 'A'
    df[ 'auth_asym_id' ] = df.label_asym_id # copy 'label_asym_id' attribute to the 'auth_asym_id' attribute.
    return df[ l_col_required ]


# ## Functions for DSSP outputs

# In[ ]:


def DSSP_Read_mkdssp_Result( dir_file ) : # 2020-07-06 17:09:09 
    ''' # 2020-12-25 18:23:24 
    Read mkdssp result file into DataFrame. Residue number in the "RESIDUE" column should be integer only.
    if directory of protein structure file is given (.pdb, .cif), run mkdssp program to read output 
    '''
    # preprocess inputs
    content = None
    if isinstance( dir_file, ( str ) ) :
        name_file = dir_file.rsplit( '/', 1 )[ 1 ]
        if '.' in name_file : # if name_file has a file extension
            str_file_extension = name_file.rsplit( '.', 1 )[ 1 ].lower( )
            if str_file_extension == 'pdb' or str_file_extension == 'cif' : # if a given file has a file extension for protein structures (an input of MKDSSP program)
                content = subprocess.run( [ 'mkdssp', '-i', dir_file ], capture_output = True ).stdout.decode( ) # run mkdssp program to read the output of mkdssp program
        if content is None : # if given file appears to not contain protein structures, read output files 
            with open( dir_file, 'r' ) as file :
                content = file.read( )
    else : # assume the given object is an opened file if type of the 'dir_file' is not string
        content = dir_file.read( )
        dir_file.close( )
        
    l_line = content.split( '\n' )
    l_l_value = list( )
    for line in l_line[ 28 : -1 ] :
        try : l_l_value.append( [ int( line[ 6 : 10 ] ), line[ 11 ], line[ 13 ], line[ 16 : 25 ], int( line[ 35 : 38 ] ), float( line[ 103 : 109 ] ), float( line[ 109 : 115 ] ) ] ) # retrieve valid line in the dssp output file
        except : continue
    df = pd.DataFrame( l_l_value, columns = [ 'residue_sequence_number', 'chain_identifier', 'amino_acid', 'structure', 'accessibility', 'phi', 'psi' ] )
    df[ 'Structure_Simple' ] = list( structure[ 0 ] if structure[ 0 ] != ' ' else 'C' for structure in df.structure ) # retrive one letter representation of a secondary structure
    
    dict_amino_acid_to_MaxASA__Tien_et_al__2013__emp__ = { 'A': 121.0, 'R': 265.0, 'N': 187.0, 'D': 187.0, 'C': 148.0, 'E': 214.0, 'Q': 214.0, 'G': 97.0, 'H': 216.0, 'I': 195.0, 'L': 191.0, 'K': 230.0, 'M': 203.0, 'F': 228.0, 'P': 154.0, 'S': 143.0, 'T': 163.0, 'W': 264.0, 'Y': 255.0, 'V': 165.0 }
    dict_amino_acid_to_maxasa = dict_amino_acid_to_MaxASA__Tien_et_al__2013__emp__ # use MaxASA__Tien_et_al__2013__emp__ data for calculating relative accessible area
    df[ 'relative_surface_accessibility' ] = list( min( 1, max( 0, record[ 'accessibility' ] / dict_amino_acid_to_maxasa[ record[ 'amino_acid' ] ] ) ) for record in df.to_dict( orient = 'record' ) ) # calculate relative surface area
    return df


# ## Functions for Designing Primers

# In[ ]:


def PRIMER_Tm( seq ) : # 2020-10-13 17:46:10 
    ''' calculate melting temperature of a given oligo sequence based on equations at http://insilico.ehu.es/tm.php?formula=basic '''
    nA, nT, nG, nC = seq.count( 'A' ), seq.count( 'T' ), seq.count( 'G' ), seq.count( 'C' )
    if len( seq ) < 14 : # if length < 14
        return (nA + nT) * 2 + (nG + nC) * 4
    else : # if length of the given oligo is > 14
        return 64.9 + 41 * (nG + nC - 16.4)/(nA + nT + nG + nC)


# ## Using Bookshelves in command line

# In[ ]:


def Bookshelves_pipeline_run_a_function_in_command_line_as_a_script( func, dir_folder, title, import_bookshelves = True ) :
    ''' # 2021-01-04 04:18:02 
    run a given function or a code snippet in command line by writing the function and current Bookshelves as a single script
    if a code snippet is given and it contains special characters like '\n', a raw string (r'some string') should be used.
    'import_bookshelves': import bookshelves scripts if True
    '''
    dir_folder_code_bookshelves = dict_setting_bookshelves[ 'codedir' ]
    if import_bookshelves : # update bookshelves python script if 'import_bookshelves' is set to True
        os.system( f"jupyter nbconvert --to script {dir_folder_code_bookshelves}Bookshelves.ipynb" )
    str_time_stamp = TIME_GET_timestamp( ) 
    dir_file_code = f'{dir_folder}{title}__Bookshelves__{str_time_stamp}.py' # directory of the code that will be written
    with open( dir_file_code, 'w' ) as newfile :
        if import_bookshelves :
            with open( f"{dir_folder_code_bookshelves}Bookshelves.py", 'r' ) as file : str_bookshelves = file.read( )
            newfile.write( str_bookshelves + '\n' )
        if callable( func ) :
            str_function = inspect.getsource( func )
            name_function = str_function.split( '(', 1 )[ 0 ].split( 'def', 1 )[ 1 ].strip( )
            newfile.write( f'\n# run a given function\n{str_function}\n# a definition of a function\n{name_function}( )\n' )
        elif isinstance( func, ( str ) ) :
            newfile.write( f'\n# run a given code snippet \n{func}\n' )
    return dir_file_code # return the name of the code file written to the storage


# ## Functions for developing web application

# In[ ]:


def Base64_Encode( dir_file_binary, dir_file_binary_base64 = None ) : # 2020-11-18 11:58:53 
    """ Perform Base64 Encoding for the given binary file.
    'dir_file_binary_base64' : default directory is 'dir_file_binary' + ".base64.txt" """
    dir_file_binary_base64 = dir_file_binary + ".base64.txt" if dir_file_binary_base64 is None else dir_file_binary_base64
    with open( dir_file_binary, 'rb' ) as file : 
        with open( dir_file_binary_base64, 'w' ) as newfile :
            newfile.write( base64.b64encode( file.read( ) ).decode( 'ascii' ) )
    return dir_file_binary_base64

# correlation 

def Correlation_Matrix_Pearsonr( arr, dtype = None ) :
    """ 
    # 2021-03-08 20:35:13 
    calculate correlation matrix of a given array (columns = cells/samples, rows = genes/observations/metrics)
    use precomputed mean-centered 'x' and 'y' divided by the matrix norm of itself and matrix multiplication to efficiently calculate correlation matrix of a gene expression matrix or other metric matrix
    """
    n = len( arr )
    if dtype is None : # default dtype np.float64
        dtype = np.float64
    # calculate mean-centered 'x' and 'y' by default
    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    arr_mean = arr.mean( axis = 1, dtype = dtype )
    arr_mean_centered = arr - arr_mean.reshape( ( len( arr_mean ), 1 ) )
    
    # calculate matrix norm of mean-centered 'x' and 'y' by default
    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    # calculate mean-centered 'x' and 'y' divided by the matrix norm of itself by default
    
    arr_mean_centered_divided_by_norm = arr_mean_centered / np.vstack( list( scipy.linalg.norm( a ) for a in arr_mean_centered ) )
    arr_r = np.matmul( arr_mean_centered_divided_by_norm, arr_mean_centered_divided_by_norm.T )

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    arr_r[ arr_r > 1 ] = 1
    arr_r[ arr_r < - 1 ] = - 1

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n / 2 - 1
    arr_prob = np.array( list( scipy.special.btdtr( ab, ab, e ) for e in ( 1 - np.abs( arr_r.ravel( ) ) ) / 2 ), dtype = dtype ).reshape( arr_r.shape )

    return arr_r, arr_prob

def Pearsonr( x, y, xm = None, ym = None, normxm = None, normym = None, xm_divided_by_normxm = None, ym_divided_by_normym = None, dtype = None, return_intermediate_results = False ) :
    """ 
    # 2021-03-08 20:35:13 
    pearsonr function of scipy.stats.pearsonr that receives precomputed values for rapid calculation of correlation values
    assumes x and y are numpy arrays
    'xm' : precomputed mean-centered x values
    'ym' : precomputed mean-centered y values
    'normxm' : precomputed matrix norm of mean-centered x values
    'normym' : precomputed matrix norm of mean-centered x values
    'xm_divided_by_normxm' : mean-centered 'x' values divided by the matrix norm of itself
    'ym_divided_by_normym' : mean-centered 'y' values divided by the matrix norm of itself
    'return_intermediate_results' : return 'xm_divided_by_normxm' and 'ym_divided_by_normym' in addition to 'r' and 'prob'
    """
    n = len(x)
    
    if dtype is None : # default dtype is dtype of the array 'x'
        dtype = x.dtype
    # calculate mean-centered 'x' and 'y' by default
    # By using `astype(dtype)`, we ensure that the intermediate calculations
    # use at least 64 bit floating point.
    if xm is None and xm_divided_by_normxm is None :
        xmean = x.mean( dtype = dtype )
        xm = x.astype( dtype ) - xmean
    if ym is None and ym_divided_by_normym is None :
        ymean = y.mean( dtype = dtype )
        ym = y.astype( dtype ) - ymean
    # calculate matrix norm of mean-centered 'x' and 'y' by default
    # Unlike np.linalg.norm or the expression sqrt((xm*xm).sum()),
    # scipy.linalg.norm(xm) does not overflow if xm is, for example,
    # [-5e210, 5e210, 3e200, -3e200]
    if normxm is None and xm_divided_by_normxm is None :
        normxm = scipy.linalg.norm( xm )
    if normym is None and ym_divided_by_normym is None :
        normym = scipy.linalg.norm( ym )
    # calculate mean-centered 'x' and 'y' divided by the matrix norm of itself by default
    if xm_divided_by_normxm is None :
        xm_divided_by_normxm = xm / normxm
    if ym_divided_by_normym is None :
        ym_divided_by_normym = ym / normym

    r = np.dot( xm_divided_by_normxm, ym_divided_by_normym )

    # Presumably, if abs(r) > 1, then it is only some small artifact of
    # floating point arithmetic.
    r = max(min(r, 1.0), -1.0)

    # As explained in the docstring, the p-value can be computed as
    #     p = 2*dist.cdf(-abs(r))
    # where dist is the beta distribution on [-1, 1] with shape parameters
    # a = b = n/2 - 1.  `special.btdtr` is the CDF for the beta distribution
    # on [0, 1].  To use it, we make the transformation  x = (r + 1)/2; the
    # shape parameters do not change.  Then -abs(r) used in `cdf(-abs(r))`
    # becomes x = (-abs(r) + 1)/2 = 0.5*(1 - abs(r)).  (r is cast to float64
    # to avoid a TypeError raised by btdtr when r is higher precision.)
    ab = n / 2 - 1
    prob = 2 * scipy.special.btdtr( ab, ab, 0.5 * (1 - abs(np.float64( r ) ) ) )

    if return_intermediate_results :
        return r, prob, xm_divided_by_normxm, ym_divided_by_normym
    else :
        return r, prob