#
#
# Hans-Rainer Kloeckner
#
# MPIfR 2022
#
#
# This is an example to get you a head start using a HDF5 files
# of the MPG telescope 
#

# -----------------------------------------------------------------------------------------
# usage:
#
# python 1_MY_FIRST_DATASET.py FILE_NAME 
#
#       the default will show data shape with timestamp in the keyword
#
# you can add as much keys to the input to check the data e.g.
#
#       python 1_MY_FIRST_DATASET.py FILE_NAME monitor weather
#
# -----------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


import h5py
import numpy as np
import sys
#
from MPG_HDF5_libs import *

#  Takes the first argument as file name
#
#
data_file = sys.argv[1]



# load the hdf5 file
#
#
obsfile = h5py.File(data_file)


# File name used 
#
print('\n\n###')
print('### data info and structure of file: ',data_file)
print('###\n\n')

# The data attributes includes all kind of information you may need for further processing
#
print('=== attributes ===\n')
obs_attrs = obsfile.attrs
#
for a in obs_attrs:
    print(' -',a,'   ',obsfile.attrs[a])
#
# ----


# The data file is group into several groups or tables that can be listed via 
#
# ['history', 'monitor', 'scan']
#   dataset     group     group
#
print('\n=== group tables ===\n')

groups_tables = obsfile.keys()

print(' -',list(groups_tables))
#
# ----


# The data file is further structured in groups which can be 
# addressed via keys similar to the dictionary mechanism in python
#
#
print('\n=== file structure and keys ===\n')

h5printstructure(obsfile)
#
# ----


# Example to optain timestamp data of the file
#
#

if len(sys.argv) == 2:
    find_key = ['timestamp']
else:
    find_key = sys.argv[2:]

print('\n\n=== extract ',find_key,'info of the file ===')

data_key = findkeys(obsfile,keys=find_key,exactmatch=False)

print('\n   === show data shape ===')

for d in data_key:
    print('\n   data[',d,'] ',obsfile[d].shape)

#
# PLEASE NOTE IF YOU WANT TO EXTRACT DATA JUST HAVE A LOOK INTO THE NEXT EXAMPLE FILE
#
 
