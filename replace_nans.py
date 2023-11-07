import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table

filename = sys.argv[1]
try:
    hdul = fits.open(filename)
except:
    print('Error reading input file')
    exit()

hdul = fits.open(filename)

table_data = hdul[1].data

print(table_data[np.where(np.isnan(table_data['DNF_Z']))]['DNF_Z'])
print(np.where(np.isnan(table_data['DNF_Z'])))
table_data['DNF_Z'][np.where(np.isnan(table_data['DNF_Z']))] = -99.0
print(table_data[np.where(np.isnan(table_data['DNF_Z']))]['DNF_Z'])
table_data['DNF_ZSIGMA'][np.where(np.isnan(table_data['DNF_ZSIGMA']))] = -99.0
table_data['DNF_ZERR_FIT'][np.where(np.isnan(table_data['DNF_ZERR_FIT']))] = -99.0

filename_out = 'out.fits'
fits.writeto(filename_out,table_data,overwrite=True)

