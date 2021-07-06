#!/usr/bin/env python
'''
This script will create a healpix density map from a catalog with a mask

Author: Nacho Sevilla
Usage: python create_hmap_from_cat.py [options]
'''

import sys, os
import numpy as np
import healpy as hp
from astropy.io import fits
from optparse import OptionParser
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("Agg")

NSIDE = 4096
RACOO = 'RA'
DECOO = 'DEC'
ZCOO = 'DNF_ZMEAN_SOF'
data_filename = ''
mask_filename = ''

def create_density_map(data_filename, mask_filename, nest):
    print('Opening',data_filename)
    hdu = fits.open(data_filename)
    print('Reading data')
    cat = hdu[1].data
    print('Opening',mask_filename)
    hdu = fits.open(mask_filename)
    print('Reading mask')
    mask = hdu[1].data 
    mask_sel = mask[mask['FRACGOOD'] > 0] # subselection of mask with fracgood > threshold
    print('Creating map')
    galmap = np.empty(hp.nside2npix(NSIDE))
    if nest:
        galmap[hp.ring2nest(NSIDE,mask_sel['HPIX'])] = 0.
    else:
        galmap[mask_sel['HPIX']] = 0.
    print('Filtering data')
    gd, = np.where((cat[ZCOO] > 0.1) & (cat[ZCOO] < 0.8)) 
    cat = cat[gd]
    th = (90-cat[DECOO])*np.pi/180.
    ph = cat[RACOO]*np.pi/180.
    print('Changing to pixels')
    pixgd = hp.ang2pix(NSIDE,th,ph,nest=True)
    print('Sorting and aggregating')
    s = np.sort(pixgd, axis=None)
    galmap_s = np.bincount(s)
    galmap[:len(galmap_s)] = galmap_s
    mask_unseen = np.full(hp.nside2npix(NSIDE),hp.UNSEEN)
    if nest:
        mask_unseen[hp.ring2nest(NSIDE,mask_sel['HPIX'])] = 1 # set to 1 those values in the selected mask to avoid them 
    else:
        mask_unseen[mask_sel['HPIX']] = 1 
    galmap[mask_unseen < 0] = hp.UNSEEN # set to UNSEEN all those values in the density map outside the mask

    return galmap

def write_density_map(galmap,outfile):
    print('Creating and writing map')
    c1 = fits.Column(name='PIXELVAL', array=galmap, format='K')
    t = fits.BinTableHDU.from_columns([c1])
    t.writeto(outfile,overwrite=True)
    print('Plotting and saving figure')
    hp.mollview(galmap,nest=True)
    plt.savefig('mollview.png')

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--infile",type="string",dest="infile",help="Input data file",default='/afs/ciemat.es/user/s/sevilla/stellarlss_project/Y3MLLOZ_weighted.fits')
    parser.add_option("--mask",type="string",dest="mask",help="Mask file",default='/afs/ciemat.es/user/s/sevilla/stellarlss_project/y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits.gz')
    parser.add_option("--mask2nest",action="store_true",dest="nest",help="Transform mask to NEST",default=False)
    parser.add_option("--outfile",type="string",dest="outfile",help="Output file",default='density_map.fits')
    (options, args) = parser.parse_args()

    data = create_density_map(options.infile,options.mask,options.nest)

    write_density_map(data,options.outfile)

    return 0
    

if __name__ == "__main__":
    main()
