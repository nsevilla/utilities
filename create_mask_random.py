#!/usr/bin/env python
'''
This script will create a random catalog on specified RA,DEC range, and optionally mask it with healpix file
Author: Nacho Sevilla
'''
import numpy as np
import healpy as hp
import random
from astropy.io import fits
from optparse import OptionParser

def create_random(npoints, ramin, ramax, decmin, decmax):

    print('Creating ',npoints,' random points')

    u1 = np.random.rand(npoints)
    u2 = np.random.rand(npoints)

    thmin = np.pi/2 - decmax*np.pi/180
    thmax = np.pi/2 - decmin*np.pi/180
    ctmin = np.cos(thmin)
    ctmax = np.cos(thmax)
    cthet = ctmin + u1*(ctmax-ctmin)

    dec = 90 - np.arccos(cthet)*180/np.pi
    ra = ramin + u2*(ramax-ramin)

    cat = [ra,dec]
    
    return cat

def mask_healpix(cat, mask_filename, nest):

    explicit = False
    hdu = fits.open(mask_filename)
    tdata = hdu[1].data
    try:
        hp.get_nside(tdata[0])
    except TypeError:
        explicit = True
        pass

    # Load mask
    good_region_mask = hp.read_map(mask_filename,partial=explicit,nest=nest)

    # Check random points against mask
    theta = (90.0 - cat[1])*np.pi/180.
    phi = cat[0]*np.pi/180.
    nside = hp.npix2nside(good_region_mask.size)
    pix = hp.ang2pix(nside,theta,phi,nest=nest)
    good, = np.where(good_region_mask[pix] == 1)

    cat_masked = np.array([cat[0][good],cat[1][good]])
    
    return cat_masked

def write_catalog(cat,output_filename):
    hdu = fits.BinTableHDU.from_columns([fits.Column(name='RA', format='D', array=cat[0]), fits.Column(name='DEC', format='D', array=cat[1])])
    hdu.writeto(output_filename,overwrite=True)

    return 0

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--output_file",type="string",dest="output_filename",help="Output random file",default='./random.fits')
    parser.add_option("--ramax",type="float",dest="ramax",help="RA maximum (degrees)",default=100)
    parser.add_option("--ramin",type="float",dest="ramin",help="RA minimum (degrees)",default=-70)
    parser.add_option("--decmax",type="float",dest="decmax",help="DEC maximum (degrees)",default=6)
    parser.add_option("--decmin",type="float",dest="decmin",help="DEC minimum (degrees)",default=-70)
    parser.add_option("--npoints",type="int",dest="npoints",help="Number of random points",default=10000000)
    parser.add_option("--use_mask",action="store_true",dest="mask_option",help="Trigger Healpix mask",default=False)
    parser.add_option("--mask_file",type="string",dest="mask_filename",help="Healpix mask",default='./mask.fits')
    (options, args) = parser.parse_args()

    if (options.decmax <= options.decmin):
        print("decmax has to be greater than decmin")
        return -1
    if (options.ramax <= options.ramin):
        print("ramax has to be greater than ramin")
        return -1
        
    random_catalog = create_random(options.npoints,options.ramin,options.ramax,options.decmin,options.decmax)

    if options.mask_option:
        print("Masking random catalog")
        try:
            random_catalog = mask_healpix(random_catalog, options.mask_filename, nest=False)
        except FileNotFoundError:
            print("Wrong file or file path for mask %s",options.mask_filename)
            return -1
   
    write_catalog(random_catalog,options.output_filename)
    return 0


if __name__ == "__main__":
    main()
