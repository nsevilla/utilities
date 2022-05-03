import healsparse as hs
import skyproj as skp
import matplotlib.pyplot as plt
from optparse import OptionParser

def main():
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--input_file",type="string",dest="infilename",help="Input file",default='map.fits')
    parser.add_option("--nside_coverage",type="int",dest="nside_coverage",help="nside coverage value",default=32)
    parser.add_option("--output_file",type="string",dest="outfilename",help="Output file",default='out.png')
    (options, args) = parser.parse_args()

    #read a map
    inmap = hs.HealSparseMap.read(options.infilename, options.nside_coverage)

    fig = plt.figure(1, figsize=(8, 6))
    fig.clf()
    ax = fig.add_subplot(111)
    m = skp.McBrydeSkyproj(ax=ax)
    _ = m.draw_hspmap(inmap)
    m.draw_inset_colorbar()
    plt.show()
    plt.savefig(options.outfilename)

if __name__ == "__main__":
    main()
