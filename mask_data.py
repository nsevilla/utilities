#!/usr/bin/env python
'''
This script will mask catalogs
Use this script to apply a healpix map to a catalog of data

Author: Nacho Sevilla
Usage: python mask_data.py
'''

import os,sys
import numpy as np
import healpy as hp
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
import fitsio
from optparse import OptionParser

nside = 4096

def mask_catalog(catalog_file,mask_file,nest):
    '''Mask input catalog'''
    print('Reading mask',mask_file)
    #mask = fitsio.read(mask_file,ext=1)['I'].ravel()
    hdu_mask = fits.open(mask_file)
    mask = hdu_mask[1].data    

    mask_exp = np.zeros(hp.nside2npix(nside))
    mask_exp[mask['HPIX']] = mask['FRACGOOD']
    if not nest:
        mask_exp = hp.reorder(mask_exp,r2n=True)
        
    print('Reading catalog',catalog_file)
    hdu = fits.open(catalog_file)
    data = hdu[1].data   

    ra = data['RA']
    dec = data['DEC']
    theta = (90.0 - dec)*np.pi/180.
    phi = ra*np.pi/180.
    pix = hp.ang2pix(nside,theta,phi,nest=True)
    print(pix)

    good = np.zeros(len(data),dtype=bool)
    for counter,p in enumerate(pix):
        if counter%100000 == 0:
            print('Masked',counter,'objects') 
    #    if ((mask_y3_a[p] > 0) and (mask_y3_b[p] == 0)):  # (mask_y1_ext[p] > 0): 
        if (mask_exp[p] > 0):
            good[counter] = True

    print('Unmasked number of objects = ',len(data))
    print('Sum of good objects = ',len(data[good]))

    return data[good]
    
def write_catalog(data,out_file):
    '''Write catalog to fits file'''
    print('Writing to file',out_file)
    hdu_masked = fits.BinTableHDU(data=data)
    hdu_masked.writeto(out_file,overwrite=True)
    return 0

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--mask",type="string",dest="mask",help="Mask file",default='/Users/nsevilla/des/data/y3y1_redmagic_comparison/y1_mask.fits.gz')
    parser.add_option("--catalog",type="string",dest="catalog",help="Catalog file",default='/Users/nsevilla/des/data/y3y1_redmagic_comparison/y3_redmagic.fits.gz')
    parser.add_option("--outfile",type="string",dest="out",help="Output file",default='masked.fits')
    parser.add_option("--nest",action="store_true",dest="nest",help="Set mask healpix coding to NEST",default=False)
    (options, args) = parser.parse_args()

    data = mask_catalog(options.catalog,options.mask,options.nest)
    write_catalog(data,options.out)
    return 0

if __name__ == "__main__":
    main()
