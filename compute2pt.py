import os,sys
import numpy as np
import treecorr
import fitsio
import matplotlib
import matplotlib.pyplot as plt
import healpy as hp
import healpix_util as hu
import time
from optparse import OptionParser

def compute_2pt(data_filename,random_filename,config,binslop,weights,csv_name):

    print('Loading data')
    start = time.time()
    if weights:
        datacat = treecorr.Catalog(data_filename, ra_col='RA', dec_col='DEC', ra_units='degrees', dec_units='degrees', w_col='weight')
    else:
        datacat = treecorr.Catalog(data_filename, ra_col='RA', dec_col='DEC', ra_units='degrees', dec_units='degrees')
    randcat = treecorr.Catalog(random_filename, ra_col='RA', dec_col='DEC', ra_units='degrees', dec_units='degrees')

    print('Starting with estimation')
    
    dd = treecorr.NNCorrelation(config,bin_slop=binslop) # var_method = 'bootstrap'
    rr = treecorr.NNCorrelation(config,bin_slop=binslop)
    dr = treecorr.NNCorrelation(config,bin_slop=binslop)

    dd.process(datacat)
    rr.process(randcat)
    dr.process(datacat,randcat)

    wtheta_counts,dwtheta_counts = dd.calculateXi(rr,dr)

    th_counts = dd.meanr
    wth_counts = wtheta_counts

    np.savetxt('theta_counts_{}.csv'.format(csv_name),th_counts)
    np.savetxt('wtheta_counts_{}.csv'.format(csv_name),wth_counts)

    end = time.time()
    print('Time for counts estimator with bin_slop {0} is {1} s'.format(binslop,end-start))

    result = [th_counts, wth_counts]

    return result

def compute_2pt_pix(datadens_filename,mask_filename,config,nside=4096):

## now the pixel estimator
    #datadens_filename = workdir+'density_map.fits'
    #mask_filename = workdir+'y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits.gz'
    print('Loading data')
    start = time.time()
    density = fitsio.read(datadens_filename,ext=1)['PIXELVAL']
    mask_pix = fitsio.read(mask_filename,ext=1)['HPIX']
    mask_fdet = fitsio.read(mask_filename,ext=1)['FRACGOOD'] 

    pixindex = np.arange(hp.nside2npix(nside))
    pixels = pixindex[density > -1] # density field already masked
    hpix = hu.HealPix('nest',nside)
    ra,dec = hpix.pix2eq(pixels)

    ngal = density[density > -1]
    ngal_avg = np.average(ngal, weights = mask_fdet)
    deltas = (ngal-ngal_avg)/ngal_avg

    print('Starting with pixel estimator')
    wtheta_pix = treecorr.KKCorrelation(config) # var_method = 'bootstrap'
    cat = treecorr.Catalog(ra=ra,dec=dec,ra_units='degrees',dec_units='degrees', k=deltas)#, w=weights)
    wtheta_pix.process(cat)
    th_pix = wtheta_pix.meanr
    wth_pix = wtheta_pix.xi

    np.savetxt('theta_pix_{}.csv'.format(nside),th_pix)
    np.savetxt('wtheta_pix_{}.csv'.format(nside),wth_pix)

    end = time.time()
    print('Time for pixel estimator with nside {0} is {1} s'.format(nside,end-time_counts))

    result = [th_pix,wth_pix]

    return result

def plot_wth(result,config,label,outfile_name,result_comp=None,label_comp=None):

    plt.xlabel('theta {}'.format(config['sep_units']))
    plt.ylabel('theta*w(theta)')
    plt.xscale('log')
    plt.yscale('linear')
    plt.xlim(config['min_sep'],config['max_sep'])
    #plt.ylim(0.0,1.0)
    plt.plot(result[0],result[0]*result[1],label=label)
    if result_comp is not None:
        plt.plot(result_comp[0],result_comp[0]*result_comp[1],label=label_comp)
    plt.legend()
    plt.savefig('thwth_{}.png'.format(outfile_name))

#    if result_comp is not None:
#        plt.figure()
#        plt.xscale('log')
#        plt.yscale('linear')
#        plt.plot(result[0],result_comp[1]-result[1])
#        plt.savefig('wth_residual_{}_{}.png'.format(label,label_comp))

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("--data_filename",type="string",dest="data_filename",help="Data file name",default="/afs/ciemat.es/user/s/sevilla/projects/stellarlss_project/data/y3maglim_bin1_std34.fits")
    parser.add_option("--random_filename",type="string",dest="random_filename",help="Random file name",default="/afs/ciemat.es/user/s/sevilla/projects/stellarlss_project/random.fits")
    parser.add_option("--min_sep",type="float",dest="min_sep",help="Minimum theta",default=2.5)
    parser.add_option("--max_sep",type="float",dest="max_sep",help="Maximum theta",default=250.0)
    parser.add_option("--nbins",type="int",dest="nbins",help="Number of bins",default=20)
    parser.add_option("--sep_units",type="string",dest="sep_units",help="Units of theta",default="arcmin")
    parser.add_option("--binslop",type="float",dest="binslop",help="bin_slop for treecorr",default=0.1) 
    parser.add_option("--weights",action="store_true",dest="weights_flag",help="Read and apply weights column",default=False) 
    parser.add_option("--plot_label",type="string",dest="plot_label",help="Label for plot in legend",default=None)
    parser.add_option("--compare",action="store_true",dest="compare_flag",help="Make comparison with another 2pcf",default=False)
    parser.add_option("--data_filename_comp",type="string",dest="data_filename_comp",help="Data file name to compare",default=None)
    parser.add_option("--random_filename_comp",type="string",dest="random_filename_comp",help="Random file name to compare",default=None)
    parser.add_option("--binslop_comp",type="float",dest="binslop_comp",help="bin_slop for treecorr for comparison test",default=0.1)
    parser.add_option("--weights_comp",action="store_true",dest="weights_comp_flag",help="Read and apply weights column to comparison data",default=False) 
    parser.add_option("--plot_label_comp",type="string",dest="plot_label_comp",help="Label for comparison plot in legend",default=None)
    parser.add_option("--outfile",type="string",dest="outfile_name",help="Label for plot",default="treecorr") 
    (options, args) = parser.parse_args()

    config = {"min_sep" : options.min_sep, 
              "max_sep" : options.max_sep, 
              "nbins" : options.nbins, 
              "sep_units" : options.sep_units, 
              "verbose" : 2}
    
    result = compute_2pt(options.data_filename,options.random_filename,config,options.binslop,options.weights_flag,"1")
    result_comp = None
    if options.compare_flag:
        result_comp = compute_2pt(options.data_filename_comp,options.random_filename_comp,config,options.binslop_comp,options.weights_comp_flag,"2")
        #result_comp = compute_2pt_pix()
    plot_wth(result,config,
             label='{}'.format(options.plot_label),
             outfile_name='{}'.format(options.outfile_name),
             result_comp=result_comp,label_comp='{}'.format(options.plot_label_comp))

if __name__ == "__main__":
    main()
