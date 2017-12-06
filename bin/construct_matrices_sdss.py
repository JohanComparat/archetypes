"""
python3.4 construct_matrices_sdss.py -0.1 0.05 X_AGN
"""
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle

version = 'v2'

# parameters of the run
zmin = float(sys.argv[1])
zmax = float(sys.argv[2])
gal_type = sys.argv[3]
# ELG, LRG, X_AGN, QSO

Nspec_max = 2000.

name = "sdss_"+gal_type+"_zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_Nlt_"+str(Nspec_max)

wlmin = 3700.
wlmax = 9000.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
d_lambda = 1. 
masterwave = n.arange(wlmin_rf, wlmax_rf, d_lambda) #0.5)

# gets the eboss elg summary file
if gal_type == "XAGN":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/specObj-dr14.fits'))[1].data
	ok_i=(
	  (cat['SOURCETYPE']=='XMMBRIGHT')|
	  (cat['SOURCETYPE']=='XMMGRIZ')|
	  (cat['SOURCETYPE']=='XMMHR')|
	  (cat['SOURCETYPE']=='XMMRED')|
	  (cat['SOURCETYPE']=='XMM_PRIME')|
	  (cat['SOURCETYPE']=='XMM_SECOND')|
	  (cat['SOURCETYPE']=='SPIDERS_PILOT')|
	  (cat['SOURCETYPE']=='SPIDERS_RASS_AGN')|
	  (cat['SOURCETYPE']=='SPIDERS_XMMSL_AGN')|
	  (cat['SOURCETYPE']=='S82X_TILE1')|
	  (cat['SOURCETYPE']=='S82X_TILE2')
	  )

if gal_type == "QSO":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/specObj-dr14.fits'))[1].data
	ok_i=(cat['ZWARNING']==0)&(cat['Z']>zmin)&(cat['Z']<=zmax)&(
	  (cat['SOURCETYPE']=='QSO')|
	  (cat['SOURCETYPE']=='QSO_EBOSS_W3_ADM')|
	  (cat['SOURCETYPE']=='QSO_GRI')|
	  (cat['SOURCETYPE']=='QSO_HIZ')|
	  (cat['SOURCETYPE']=='QSO_RIZ')|
	  (cat['SOURCETYPE']=='QSO_VAR')|
	  (cat['SOURCETYPE']=='QSO_VAR_FPG')|
	  (cat['SOURCETYPE']=='QSO_VAR_LF')|
	  (cat['SOURCETYPE']=='QSO_VAR_SDSS')|
	  (cat['SOURCETYPE']=='QSO_WISE_FULL_SKY')|
	  (cat['SOURCETYPE']=='QSO_WISE_SUPP')|
	  (cat['SOURCETYPE']=='QSO_XD_KDE_PAIR')|
	  (cat['SOURCETYPE']=='RADIO_2LOBE_QSO')|
	  (cat['SOURCETYPE']=='WISE_BOSS_QSO')|
	  (cat['SOURCETYPE']=='TDSS_FES_HYPQSO' ) )
	
if gal_type == "ELG":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/ELG/eboss21/cats/eboss21.v5_10_7.latest.fits'))[1].data
	zQ_min = 2.
	ok_i=(cat['Z_ZQ']>=zQ_min)&(cat['Z']>zmin)&(cat['Z']<=zmax)

if gal_type == "LRG":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/specObj-dr14.fits'))[1].data
	ok_i=(cat['Z']>zmin)&(cat['Z']<=zmax)&((cat['SOURCETYPE']=='LRG')|(cat['SOURCETYPE']=='FAINT_HIZ_LRG')|(cat['SOURCETYPE']=='HIZ_LRG')|(cat['SOURCETYPE']=='LRG_ROUND3'))


# mask for the sky lines
maskLambda = n.loadtxt(os.path.join(os.environ['GIT_ARCHETYPES'],'data',"dr12-sky-mask.txt"), unpack=True)

get_path_to_spectrum_v5_10_0 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_0', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )
get_path_to_spectrum_v5_10_7 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_7', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )
get_path_to_spectrum_26 = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', '26', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )


def get_spec(path_to_spectrum, maskLambda=maskLambda):
	"""
	for a spectrum, returns plate, mjd, fiberid
	"""
	hdulist = fits.open(path_to_spectrum)
	wavelength = 10**hdulist[1].data['loglam']
	flux = hdulist[1].data['flux']
	ivar = hdulist[1].data['ivar']
	ratio = n.min(abs(10000.*n.log10(n.outer(wavelength, 1./maskLambda))), axis=1)
	margin = 1.5
	veto_sky = ratio <= margin
	selection = (veto_sky) & (ivar<=0) & (flux<0.)& (n.isinf(ivar)) & (n.isinf(flux))
	flux[selection] = n.zeros_like(ivar[selection])
	ivar[selection] = n.zeros_like(ivar[selection])
	return wavelength, flux, ivar #, selection


# subselection
nGal_all = len(ok_i.nonzero()[0])
print("total", nGal_all)
if nGal_all > Nspec_max:
    rds = n.random.rand(len(cat))
    ok = (ok_i)&(rds < Nspec_max / nGal_all)
else :
    ok = ok_i
    
nGal = len(ok.nonzero()[0])
#z=cat['Z'][ok]


print("run name:",name)
print("considering", nGal)

allflux = n.zeros( ( nGal, len(masterwave)) )
allivar = n.zeros( ( nGal, len(masterwave)) )

print("gets the spectra to construct the matrix")
for index in range(nGal) :# = 0
	path_1 = get_path_to_spectrum_v5_10_0(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	path_2 = get_path_to_spectrum_26(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	path_3 = get_path_to_spectrum_v5_10_7(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	choice = n.array([os.path.isfile(path_1), os.path.isfile(path_2), os.path.isfile(path_3)])
	pathes = n.array([path_1, path_2, path_3])
	print(pathes, choice)
	if choice.any():
		path = pathes[choice][0]
		print( "exists" )
		redshift = cat['Z'][ok][index]
		wl, fl, iv = get_spec(path)
		#wl_rf = wl / (1+redshift)
		w_new = masterwave*(1+redshift)
		itp_fl = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,wl*fl,0.,0.)) )
		itp_iv = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,iv,0.,0.)) )
		try :
			f_new = itp_fl(w_new)/w_new
			iv_new = itp_iv(w_new)
			allflux[index] = f_new
			allivar[index] = iv_new
		except(ValueError):	
			pass
		
out_dir = os.path.join( os.environ['OBS_REPO'], 'archetypes', version, gal_type )
if os.path.isdir(out_dir)==False:
  os.system('mkdir -p '+out_dir)
  
n.savetxt(os.path.join(out_dir, "allflux_"+name+".txt")   , allflux)
n.savetxt(os.path.join(out_dir, "allivar_"+name+".txt")   , allivar)
n.savetxt(os.path.join(out_dir, "masterwave_"+name+".txt"), masterwave)


class ObjIds:
    def __init__(self, plate, mjd, fiberid):
        self.plate = plate
        self.mjd = mjd
        self.fiberid = fiberid
        
f = open(os.path.join(out_dir, "specObj_"+name+".pkl"), 'wb')
obj = ObjIds(cat['PLATE'][ok], cat['MJD'][ok], cat['FIBERID'][ok])
pickle.dump(obj, f)
f.close()
