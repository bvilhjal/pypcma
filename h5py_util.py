"""
A couple of functions to convert dictionaries to HD5 files.

Use with cauthion
"""

import h5py

def dict_to_hdf5(input_dict, hdf5_group):
    """
    Recursively constructs HDF5 file groups from dictionaries.  
    
    Assumes that the data can be undestood by h5py create_dataset function.
    """
    for key in input_dict:
        if isinstance(input_dict[key], dict):
            new_hdf5_group = hdf5_group.create_group(key)
            dict_to_hdf5(input_dict[key], new_hdf5_group)
        else:
            hdf5_group.create_dataset(key, data=input_dict[key])
            
def hdf5_to_dict(hdf5_group):
    output_dict = {}
    for key in hdf5_group.keys():
        if isinstance(hdf5_group[key], h5py._hl.group.Group):
            output_dict[key] = hdf5_to_dict(hdf5_group[key])
        else:
            output_dict[key] = hdf5_group[key][...]
    return output_dict
    