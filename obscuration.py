#!/usr/bin/env python
'''
This script will plot the amount of obscuration in a MAIN sample from a BRIGHT sample
Author: Nacho Sevilla, based on work from Jelena Aleksic
Usage: python obscuration.py --bright_catalog [CATALOG WITH OBSCURING OBJECTS] --main_catalog [MAIN CATALOG] 
'''
import numpy as np
import os,sys
import sklearn
from sklearn.neighbors import NearestNeighbors as NN
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("Agg")
#matplotlib.style.use('des_dr1')
from descolors import BAND_COLORS 
from astropy.io import fits
from astropy.io import ascii
from optparse import OptionParser

RACOO1 = 'ra_deg_cont'
DECOO1 = 'dec_deg_cont'
#RACOO1 = 'Ra'
#DECOO1 = 'Dec'
#RACOO2 = 'ra_deg_cont'
#DECOO2 = 'dec_deg_cont'
#RACOO2 = 'Ra'
#DECOO2 = 'Dec'
RACOO2 = 'RA'
DECOO2 = 'DEC'

def ring(i,intr,extr):
    area = np.pi*(extr*extr-intr*intr)
    return area

def compute_obscuration(filename_bright, filename_main, distance, nbins, cutname):

    bins,step = np.linspace(0., distance, nbins+1, retstep = True)
    midbins = bins+0.5
    midbins = midbins[0:-1]

    #real data
    hdulist = fits.open(filename_bright,memmap=True)
    tdata_bright = hdulist[1].data
    hdulist = fits.open(filename_main,memmap=True)
    tdata_main = hdulist[1].data

    # Additional cut to be applied to the bright sample
    cut     = [[50, 100], [100, 200], [200, 500], [500, 1000]]
    #cut     = [[0.0005, 0.001], [0.001, 0.002]]
    ncut = len(cut)

    # BRIGHT
    #cutB = (tdata_bright['flag_i1'] < 1) # dummy cut
    cutB = (tdata_bright[RACOO1] > 310) & (tdata_bright[RACOO1] < 330) & (tdata_bright[DECOO1] > -61) & (tdata_bright[DECOO1] < -50) # dummy cut
    xB = tdata_bright[RACOO1][cutB]*60. #degrees to arcmin
    yB = tdata_bright[DECOO1][cutB]*60. #degrees to arcmin
    #star = np.vstack((xS,yS)).T
    var = tdata_bright[cutname][cutB]
    nB = len(tdata_bright[cutB])

    # MAIN
    #cutM = (tdata_bright['flag_i1'] < 1) # dummy cut
    #cutM = (tdata_main['Included'] > -1) # dummy_cut
    cutM =   (tdata_main[RACOO2] > -10000) #dummy_cut
    xM = tdata_main[RACOO2][cutM]*60. #degrees to arcmin
    yM = tdata_main[DECOO2][cutM]*60. #degrees to arcmin
    nM = len(tdata_main[cutM])

    print(' * Applying additional cut to bright sample in', cutname)
    xBbin, yBbin, nBbin = [], [], []
    for i in range(ncut):
        xBbin.append(xB[(cut[i][0] <= var) & (var < cut[i][1])])
        yBbin.append(yB[(cut[i][0] <= var) & (var < cut[i][1])])
        nBbin.append(len(xBbin[i]))
        print('     ** Cut', i+1, ': (', cut[i][0], '<=', cutname, '<', cut[i][1], ') --> Bright:', nBbin[i])

    bright_sample, main_sample = [], []
    for i in range(ncut):
        tmpB = []
        for j in range(nBbin[i]): tmpB.append([xBbin[i][j], yBbin[i][j]])
        bright_sample.append(np.asarray(tmpB))
        del tmpB
    for i in range(nM):
        main_sample.append([xM[i], yM[i]])
    main_sample = np.asarray(main_sample)

    print('Calculating NN')
    neigh = NN(radius = distance+1, metric = 'euclidean')
    neigh.fit(main_sample)
    distB = [] #distance from central bright object
    for i in range(ncut):
        dist, ind = neigh.radius_neighbors(bright_sample[i])
        distB.append(dist)

    print('Calculating histos')
    hist,histerr = [],[]
    area = np.zeros((nbins))

    colors = [BAND_COLORS['u'],BAND_COLORS['g'],BAND_COLORS['r'],BAND_COLORS['i'],BAND_COLORS['z'],BAND_COLORS['Y']] #these are colors defined for DES plots, to be used in general
    plt.figure()
    for i in range(ncut):
        d, derr, cnt = np.zeros((nbins)), np.zeros((nbins)), np.zeros((nbins))
        #print(nBbin[i],nbins)
        for j in range(nBbin[i]):
            tmphist, _ = np.histogram(distB[i][j], bins)
            d += tmphist
            derr += tmphist
            cnt += 1
        for k in range(nbins):
            area[k] = ring(k+1,_[k],_[k+1]) # i=0 will give you funny results as area() is defined
        d = d/cnt/area
        derr = np.sqrt(derr)/cnt/area
        hist.append(d/d[nbins-1])
        histerr.append(derr/d[nbins-1])
        print(d)
        print(midbins)
        print(hist)
        print(derr)
        print(len(midbins),len(hist))
        plt.errorbar(midbins, hist[i], yerr=histerr[i], fmt='o', label=cutname+str(cut[i]), color=colors[i])
    plt.ylabel('Relative abundance of sources')
    plt.xlabel('Distance from bright source (arcmin)') 
    plt.legend()
    plt.savefig('obscuration.png')

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--bright_catalog",type="string",dest="filename_bright",help="Bright catalog file",default='/Users/nsevilla/emu/data/pilot_islands.fits')
    parser.add_option("--main_catalog",type="string",dest="filename_main",help="Main catalog file",default='/Users/nsevilla/emu/data/pilot_islands.fits')
    parser.add_option("--distance",type="float",dest="distance",help="Maximum distance in arcmin",default=40)
    parser.add_option("--nbins",type="int",dest="nbins",help="Number of bins",default=10)
    parser.add_option("--cutname",type="string",dest="cutname",help="Variable in which to bin the bright catalog",default='flux_int')
    (options, args) = parser.parse_args()

    compute_obscuration(options.filename_bright, options.filename_main, options.distance, options.nbins, options.cutname)
    return 0

if __name__ == "__main__":
    main()
