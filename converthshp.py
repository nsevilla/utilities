import healsparse as hs
import healpy as hp
from optparse import OptionParser

def main():
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--input_file",type="string",dest="infilename",help="Input file",default='map.fits')
    parser.add_option("--nside_coverage",type="int",dest="nside_coverage",help="nside coverage value",default=32)
    parser.add_option("--output_healsparse",action="store_true",dest="ouths",help="Toggle healsparse output",default=False)
    parser.add_option("--output_file",type="string",dest="outfilename",help="Output file",default='out.fits')
    (options, args) = parser.parse_args()
    
    #read a map
    inmap = hs.HealSparseMap.read(options.infilename, options.nside_coverage)

    #write maps
    if options.ouths:
        inmap.write(options.outfilename, clobber=True)
    else:
        #hpmap_converted = inmap[:]
        hpmap_converted = inmap.generate_healpix_map(nside=4096, reduction='mean')
        hp.write_map(options.outfilename, hpmap_converted, overwrite=True)

if __name__ == "__main__":
    main()
