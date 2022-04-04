# load internal module
from biobookshelf.main import *

def shared_array_to_numpy_array( mp_arr, dtype ) :
    ''' return a numpy representation of a shared array so that numpy operation can be performed '''
    return np.frombuffer( mp_arr.get_obj( ), dtype = dtype )