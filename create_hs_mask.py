import healsparse as hs
import numpy as np
import pandas as pd

def box_mask(ra_corners, dec_corners):

    poly = hs.Polygon(ra=ra_corners, dec=dec_corners, value=1)
    smap_poly = poly.get_map(nside_coverage=32, nside_sparse=4096, dtype=np.int16)
    return smap_poly

def main():

    filename = '../radec_corners.csv'
    tiles = pd.read_csv(filename,header=None) 
    nb_corners  = int(len(tiles.columns)/2) #tiles.columns should be even
    mask = hs.HealSparseMap.make_empty(nside_coverage=32, nside_sparse=4096, dtype=np.int16)

    for i in np.arange(len(tiles)):   
        ra_corners = np.array(tiles.iloc[i,:nb_corners])
        dec_corners = np.array(tiles.iloc[i,nb_corners:nb_corners+4])
        single_mask = box_mask(ra_corners, dec_corners)
        mask = hs.or_union([mask,single_mask]) 
    
    mask.write('bad_regions_gold20.hs', clobber=False)
  
    return 0

if __name__ == "__main__":
    main()
