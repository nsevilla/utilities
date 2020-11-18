import os,sys
import numpy as np
import healpy as hp
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
import fitsio

nside = 4096

## use this script to apply a healpix map to a catalog of data

print 'Reading mask...'
#mask_y1_ring = fitsio.read('/pool/pcae75/data1/des/y1a1/gold/y1a1_gold_masks/y1a1_wide_combined_mask_bad3_depth22_fracdet0.8_4096.fits',ext=1)['I'].ravel()
#mask_y1 = hp.reorder(mask_y1_ring, r2n=True)
mask_y3_a = fitsio.read('/pool/pcae75/data1/des/y3_validation/gold_v2.2/y3a2_footprint_griz_1exp_v2.0.fits.gz',ext=1)['I'].ravel()
mask_y3_b = fitsio.read('/pool/pcae75/data1/des/y3_validation/gold_v2.2/y3a2_foreground_mask_v2.1.fits.gz',ext=1)['I'].ravel()
#### obsolete
#hdu = fits.open('/data1/des/y3_validation/gold_v2.2/y3a2_footprint_griz_1exp_v2.0.fits.gz')
#mask_y3 = hdu[1].data
#hdu = fits.open('/pool/pcae75/data1/des/y1a1/lss_masks/mask_1.0.2_fdetgriz08_badle3_i22_4096_ring.fits')
#mask_y1 = hdu[1].data
#mask_y1 = hp.read_map('/pool/pcae75/data1/des/y1a1/gold/y1a1_gold_masks/y1a1_wide_combined_mask_bad3_depth22_fracdet0.8_4096.fits')
#mask_y1 = hp.read_map('/data1/des/y1a1/lss_masks/mask_1.0.2_fdetgriz08_badle3_i22_4096_ring.fits')

print 'Reading data...'
#hdu = fits.open('/pc/desdsk01/des/sg_challenge/external_catalogs/gaia_dr2_spt_2669.fits')
#hdu = fits.open('/data1/des/y1a1/gold/208_adagold_adamof_galgold_sgvalidation.fits',memmap=True)
#hdu = fits.open('/data1/des/y1a1/gold/lss_baored_sgpsep_0610_unmasked.fits',memmap=True)
#hdu = fits.open('/pool/pcae75/data1/des/y3_validation/y3y1_density_tests/y1gold_galaxysel_s82_ragt315_nogoldmodest.fits')
#hdu = fits.open('/pool/pcae75/data1/des/y3_validation/snx3_hscw02_overlap.fits')
#hdu = fits.open('/pool/pcae75/data1/des/y3_validation/deep_fields/snx3_df_masked.fits')
hdu = fits.open('/pool/pcae75/data1/des/y3_validation/hsc/hsc_w02.fits')
#hdu = fits.open('/pool/pcae75/data1/des/y3_validation/y1gold_hscw05_completeness_v2.fits')
data = hdu[1].data
#catalogfile = '/pool/pcae75/data3/des/desemu_project/Continuum_Island_Catalogue.csv'
#data = Table.read(catalogfile, format='ascii.csv')

#ra = data['ra']
#dec = data['dec']
#ra = data['ALPHAWIN_J2000']
#dec = data['DELTAWIN_J2000']
ra = data['RA']
dec = data['DEC']
theta = (90.0 - dec)*np.pi/180.
phi = ra*np.pi/180.
pix = hp.ang2pix(nside,theta,phi,nest=True)

#mask_y1_ext = np.zeros(12*nside*nside)
#print 'Filling mask_y1_ext',len(mask_y1)
#for i in range(len(mask_y1)):
#    if i % 10000 == 0:
#	print i
#    mask_y1_ext[mask_y1['pixel_id'][i]] = mask_y1['fraction'][i] 

good = np.zeros(len(data),dtype=bool)
for counter,p in enumerate(pix):
    if counter%100000 == 0:
        print 'Masked',counter,'objects' 
    if ((mask_y3_a[p] > 0) and (mask_y3_b[p] == 0)):  # (mask_y1_ext[p] > 0): 
#    if (mask_y1[p] > 0):
#    if (mask_y1_ext[p] > 0):
        good[counter] = True

print 'Unmasked number of objects = ',len(data)
print 'Sum of good objects = ',len(data[good])#,sum(good)
print 'Writing to file...'

#data[good].write('/pool/pcae75/data3/des/desemu_project/Continuum_Island_Catalogue_DESmasked.csv', format='csv')
hdu_masked = fits.BinTableHDU(data=data[good])
hdu_masked.writeto('/pool/pcae75/data1/des/y3_validation/w02_hsc_y3goldmasked.fits',clobber=True)



