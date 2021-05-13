import sys, os
import numpy as np
import healpy as hp
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("Agg")
from astropy.io import fits

NSIDE = 128
RACOO = 'ra'
DECOO = 'dec'
ZCOO = 'zredmagic'
filename = '/Users/nsevilla/des/data/red.fits'

def main():
    print('Opening',filename)
    hdu = fits.open(filename)
    print('Reading data')
    cat = hdu[1].data
    print('Filtering data')
    gd, = np.where((cat[ZCOO] > 0.1) & (cat[ZCOO] < 0.8))
    cat = cat[gd]
    th = (90-cat[DECOO])*np.pi/180.
    ph = cat[RACOO]*np.pi/180.
    print('Changing to pixels')
    pixgd = hp.ang2pix(NSIDE,th,ph,nest=True)
    galmap = np.zeros(hp.nside2npix(NSIDE))
    print('Sorting and aggregating')
    s = np.sort(pixgd, axis=None)
    galmap_s = np.bincount(s)
    galmap[:len(galmap_s)] = galmap_s
    print('Creating and writing map')
    c1 = fits.Column(name='PIXELVAL', array=galmap, format='K')
    t = fits.BinTableHDU.from_columns([c1])
    t.writeto('map.fits')
    print('Plotting and saving figure')
    hp.mollview(galmap,nest=True)
    plt.savefig('mollview.png')

if __name__ == "__main__":
    main()
