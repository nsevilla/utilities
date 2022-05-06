import healsparse as hs
import healpy as hp
import skyproj as skp
import matplotlib.pyplot as plt
from optparse import OptionParser

def main():
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--input_file",type="string",dest="infilename",help="Input file",default='y6_MAPNAME_STATNAME_BANDNAME_4096_NEST.fits')
    parser.add_option("--nside_coverage",type="int",dest="nside_coverage",help="nside coverage value",default=32)
    parser.add_option("--nest",action="store_true",dest="isnest",help="Toggle NEST to True",default=False)
    parser.add_option("--healpix",action="store_true",dest="ishealpix",help="Toggle healpix format to True",default=False)
    parser.add_option("--output_file",type="string",dest="outfilename",help="Output file",default='out.png')
    (options, args) = parser.parse_args()

    #read a map
    print('Reading',options.infilename)
    if options.ishealpix:
        inmap = hp.read_map(options.infilename, nest = options.isnest)
    else:
        inmap = hs.HealSparseMap.read(options.infilename, options.nside_coverage)

    fig = plt.figure(1, figsize=(8, 6))
    fig.clf()
    ax = fig.add_subplot(111)
    m = skp.McBrydeSkyproj(ax=ax)

    if options.ishealpix:
        _ = m.draw_hpxmap(inmap,nest=options.isnest)
    else:        
        _ = m.draw_hspmap(inmap)
    m.draw_inset_colorbar()
    #plt.show()
    plt.savefig(options.outfilename)

if __name__ == "__main__":
    main()
