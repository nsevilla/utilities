#!/usr/bin/env python
'''
This script will mask catalogs 
Use this script to apply a healpix map to a catalog 

Author: Nacho Sevilla
Usage: python mask_data.py [options]
'''

import os,sys
import numpy as np
import healpy as hp
import healsparse as hsp
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
import fitsio
from optparse import OptionParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('Agg')

nside = 4096

def mask_data(data_file,mask_file,nest,healsparse):
    '''Mask input catalog'''
    print('Reading mask',mask_file)

    if healsparse:
        mask_exp = hsp.HealSparseMap.read(mask_file)
    else:
        #mask = fitsio.read(mask_file,ext=1)['I'].ravel()
        hdu_mask = fits.open(mask_file)
        mask = hdu_mask[1].data    
        mask_exp = np.full(hp.nside2npix(nside),hp.UNSEEN)
        #mask_exp[mask['HPIX']] = mask['FRACGOOD']
        mask_exp[mask['PIXEL']] = mask['I']
        if not nest:
            mask_exp = hp.reorder(mask_exp,r2n=True)
        
    print('Reading catalog',data_file)
    hdu = fits.open(data_file)
    data = hdu[1].data   

    #print(np.min(mask_exp))

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
        if (mask_exp[p] > 0):
            good[counter] = True

    print('Unmasked number of objects = ',len(data))
    print('Sum of good objects = ',len(data[good]))

    return data[good]
    
def write_data(data,out_file):
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
    parser.add_option("--infile",type="string",dest="infile",help="Input data file",default='/Users/nsevilla/des/data/y3y1_redmagic_comparison/y3_redmagic.fits.gz')
    parser.add_option("--outfile",type="string",dest="outfile",help="Output file",default='masked.fits')
    parser.add_option("--nest",action="store_true",dest="nest",help="Set mask healpix coding to NEST",default=False)
    parser.add_option("--healsparse",action="store_true",dest="healsparse",help="Set input mask as in healsparse format",default=False)
    (options, args) = parser.parse_args()

    data = mask_data(options.infile,options.mask,options.nest,options.healsparse)

    write_data(data,options.outfile)

    return 0

if __name__ == "__main__":
    main()
