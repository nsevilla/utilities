import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits

def selector(x_data,y_data,x_array,y_array):
    
    f_array = interp1d(x_array,y_array)
    selection = (y_data < f_array(x_data))
    selected = float(len(selection[selection == True]))
    total = float(len(x_data))
    print(selected, total)
    return selected/total

def load_data(filename,xname,yname,zname,zmin,zmax):

     hdu = fits.open(filename)
     filedata = hdu[1].data
     sel = (filedata[zname]>zmin) & (filedata[zname]<zmax)    
     return (filedata[xname][sel],filedata[yname][sel])
    
def get_threshold(tablename):

    if tablename == "gi_zKs_sgsep":
        print("Loading definition",tablename)
        #xarr = [-1.5, 0.5, 1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 8.0]
        #yarr = [0.6, 0.6, 0.75, 0.9, 1.1, 1.1, 1.2, 1.15, 1.2, 1.2]
        xarr = [-1.5,0.5,1.25,1.5,2.,8.]
        yarr = [0.6,0.6,0.85,0.95,1.15,1.15]
    else:
        print("Separator definition",tablename,"not found") 
        return -1        

    return (xarr, yarr)

def main():

    filename = "/scratch/sevilla/maglim_for_selector.fits"
    (xarr, yarr) = get_threshold("gi_zKs_sgsep")
    xname = "gi"
    yname = "J-Kspnt"     
    zbins = [0.2, 0.4, 0.55, 0.7, 0.85, 0.95, 1.05]
    for z in range(6):
        data = load_data(filename,xname,yname,"Z_MEAN",zbins[z],zbins[z+1])
        contamination = selector(data[0],data[1],xarr,yarr)
        print("Star contamination in bin",zbins[z],"-",zbins[z+1],"is",contamination*100,"%")    

if __name__ == "__main__":
    main()
