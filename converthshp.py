import healsparse as hs
import healpy as hp
from optparse import OptionParser

def main():
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--input_file",type="string",dest="infilename",help="Input file",default='/pool/cosmo01_data1/des/y6_sp_maps/mangle_maps/SPmaps/band_g/y6a2_g_o.4096_t.32768_AIRMASS.MAX_EQU.fits.fz')
    parser.add_option("--nside_coverage",type="int",dest="nside_coverage",help="nside coverage value",default=32)
    parser.add_option("--nside_out",type="int",dest="nside_out",help="nside of output file",default=4096)
    parser.add_option("--nest",action="store_true",dest="isnest",help="Toggle NEST to True",default=False)
    parser.add_option("--healpix",action="store_true",dest="ishealpix",help="Toggle healpix format to True",default=False)
    parser.add_option("--output_healsparse",action="store_true",dest="ouths",help="Toggle healsparse output",default=False)
    parser.add_option("--output_dir",type="string",dest="outdir",help="Output directory",default='./')
    parser.add_option("--output_file_mapname",type="string",dest="mapname",help="Output file map name",default='MAPNAME')
    parser.add_option("--output_file_statname",type="string",dest="statname",help="Output file stat name",default='STATNAME')
    parser.add_option("--output_file_bandname",type="string",dest="bandname",help="Output file band name",default=None)
    (options, args) = parser.parse_args()
    
    print('Selected NEST',options.isnest)
    print('nside',options.nside_out)

    #read a map
    print('Reading map',options.infilename)
    print('Selected input format as Healpix',options.ishealpix)
    if options.ishealpix:
        inmap = hp.read_map(options.infilename, nest = options.isnest)
    else:
        inmap = hs.HealSparseMap.read(options.infilename, options.nside_coverage)

    #build string with file name
    if options.isnest:
        nestring = 'NEST'
    else:
        nestring = 'RING' 

    if options.bandname is None:
        outfilename_noext = '_'.join(('y6',options.mapname,options.statname,str(options.nside_out),nestring))
    else:
        outfilename_noext = '_'.join(('y6',options.mapname,options.statname,options.bandname,str(options.nside_out),nestring))

    if options.ouths:
        extension = '.hs'
    else:
        extension = '.fits'

    outfilename = options.outdir + outfilename_noext + extension
 
    #write maps
    print('Writing map',outfilename)
    if options.ishealpix: #input is healpix
        if options.ouths: #output in healsparse format
            conv_map = hs.HealSparseMap(nside_coverage=options.nside_coverage, healpix_map=inmap)
            conv_map.write(outfilename, clobber=True)  
        else: #output in healpix
            if options.isnest:
                order = 'NESTED'
            else:
                order = 'RING'
            lores_map = hp.ud_grade(inmap, options.nside_out, order_in = order, order_out = order)
            hp.write_map(outfilename, lores_map, nest = options.isnest, overwrite=True)
    else: #input is healsparse
        if options.ouths: #output in healsparse format
            inmap.write(outfilename, clobber=True)    
        else: #output in healpix
            conv_map = inmap.generate_healpix_map(nside=options.nside_out, reduction='mean', nest = options.isnest)
            hp.write_map(outfilename, conv_map, nest = options.isnest, overwrite=True)

if __name__ == "__main__":
    main()
