import sys
import numpy as np
from astropy.io import fits

filename = sys.argv[1]
try:
    hdul = fits.open(filename)
except:
    print('Error reading input file')
    exit()

table_data = hdul[1].data
for colname in hdul[1].columns.names:
    data = table_data[colname]
    has_nan = np.isnan(data).any()
    if(has_nan):
        print(colname,has_nan)
        nan_indices = np.argwhere(np.isnan(data))
        print('Number of NaNs:',len(nan_indices))
        print(nan_indices)
        print(table_data[nan_indices[0]]) 
