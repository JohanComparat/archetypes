"""
python3.4 construct_matrices_sdss.py -0.1 0.05 X_AGN
"""
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle
import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p

version = 'v5'

# parameters of the run
zmin = 1. # float(sys.argv[1])
zmax = 2. # float(sys.argv[2])
gal_type = 'cosmosagn' # sys.argv[3]
# ELG, LRG, X_AGN, QSO

Nspec_max = 3000

name = "deimos_"+gal_type+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))

SN_threshold = 0.1#0.5
wlmin = 6000.
wlmax = 9500.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
d_lambda = 0.2 
masterwave = n.arange(wlmin_rf, wlmax_rf, d_lambda) #0.5)

cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'COSMOS', 'DEIMOS_catalog_27MAR2018.fits'))[1].data

ok=(cat['z']>zmin)&(cat['z']<=zmax)
print(len(cat['z'][ok]))
print(n.histogram(cat['Z'][ok],bins=n.arange(0,3.1,0.5)))
nGal = len(ok.nonzero()[0])

#cat['Mask']
#cat['Slt']
#cat['Name']

get_path_to_spectrum = lambda slit, name : glob.glob(os.path.join( os.environ['OBS_REPO'], 'COSMOS', 'spec1d.*-*.' +str(slit).zfill(3) +"."+name+"_acal.fits" ))[0]


def get_spec(path_to_spectrum):
	hdulist = fits.open(path_to_spectrum)
	wavelength = hdulist[1].data['LAMBDA'][0]
	flux = hdulist[1].data['FLUX'][0]
	ivar = hdulist[1].data['IVAR'][0]
	sel = (flux>0)&(ivar>0)
	return wavelength[sel], flux[sel], ivar[sel] #, selection

allflux = n.zeros( ( nGal, len(masterwave)) )
allivar = n.zeros( ( nGal, len(masterwave)) )

print("gets the spectra to construct the matrix")
for index in range(nGal) :# = 0
	try :
		print(index, cat['Mask'][ok][index], cat['Slt'][ok][index], cat['Name'][ok][index])
		path = get_path_to_spectrum(cat['Slt'][ok][index], cat['Name'][ok][index])
		print( path )
		redshift = cat['z'][ok][index]
		wl, fl, iv = get_spec(path)
		p.figure(0,(10,5))
		p.plot(wl,fl)
		p.savefig(path+".png")
		p.clf()
		w_new = masterwave*(1+redshift)
		itp_fl = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,80000.)), n.hstack((0.,0.,wl*fl,0.,0.)) )
		itp_iv = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,80000.)), n.hstack((0.,0.,iv,0.,0.)) )
		f_new = itp_fl(w_new)/w_new
		iv_new = itp_iv(w_new)
		allflux[index] = f_new
		allivar[index] = iv_new
	except(ValueError, IndexError):	
		print('error', index, cat['Slt'][ok][index], cat['Name'][ok][index])
		
out_dir = os.path.join( os.environ['OBS_REPO'], 'archetypes', version, gal_type )
if os.path.isdir(out_dir)==False:
  os.system('mkdir -p '+out_dir)
  
n.savetxt(os.path.join(out_dir, "allflux_"+name+".txt")   , allflux)
n.savetxt(os.path.join(out_dir, "allivar_"+name+".txt")   , allivar)
n.savetxt(os.path.join(out_dir, "masterwave_"+name+".txt"), masterwave)


#class ObjIds:
    #def __init__(self, plate, mjd, fiberid):
        #self.plate = plate
        #self.mjd = mjd
        #self.fiberid = fiberid
        
#f = open(os.path.join(out_dir, "specObj_"+name+".pkl"), 'wb')
#obj = ObjIds(cat['PLATE'][ok], cat['MJD'][ok], cat['FIBERID'][ok])
#pickle.dump(obj, f)
#f.close()
