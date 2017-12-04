import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits

# gets the eboss elg summary file
path_2_eboss21_cata = os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/specObj-dr14.fits')
cat = fits.open(path_2_eboss21_cata)[1].data

# mask for the sky liens
maskLambda = n.loadtxt(os.path.join(os.environ['GIT_ARCHETYPES'],'data',"dr12-sky-mask.txt"), unpack=True)

get_path_to_spectrum = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_0', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )

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


# parameters of the run
zmin = float(sys.argv[1])
zmax = float(sys.argv[2])
wlmin = 3600.
wlmax = 9900.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
masterwave = n.arange(wlmin_rf, wlmax_rf, 0.5)

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
  (cat['SOURCETYPE']=='XMMBRIGHT')|
  (cat['SOURCETYPE']=='XMMGRIZ')|
  (cat['SOURCETYPE']=='XMMHR')|
  (cat['SOURCETYPE']=='XMMRED')|
  (cat['SOURCETYPE']=='XMM_PRIME')|
  (cat['SOURCETYPE']=='XMM_SECOND')|
  (cat['SOURCETYPE']=='SPIDERS_PILOT')|
  (cat['SOURCETYPE']=='SPIDERS_RASS_AGN')|
  (cat['SOURCETYPE']=='SPIDERS_XMMSL_AGN')|
  (cat['SOURCETYPE']=='TDSS_FES_HYPQSO')
  )

nGal_all = len(ok_i.nonzero()[0])
print("total", nGal_all)
if nGal_all > 5000.:
    rds = n.random.rand(len(cat))
    ok = (ok_i)&(rds < 5000. / nGal_all)
else :
    ok = ok_i
    
nGal = len(ok.nonzero()[0])
galaxytype=cat['CLASS'][ok]
z=cat['Z'][ok]
zgood = n.ones_like(z)


name = "ebossQSO_zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_Ngal_"+str(nGal)
print("run name:",name)
print("considering", nGal)

allflux = n.zeros( ( nGal, len(masterwave)) )
allivar = n.zeros( ( nGal, len(masterwave)) )

print("gets the spectra to construct the matrix")
for index in range(nGal) :# = 0
	path = get_path_to_spectrum(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	print(path)
	if os.path.isfile(path):
		print( "exists" )
		redshift = cat['Z'][ok][index]
		wl, fl, iv = get_spec(path)
		#wl_rf = wl / (1+redshift)
		w_new = masterwave*(1+redshift)
		itp_fl = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,wl*fl,0.,0.)) )
		itp_iv = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,iv,0.,0.)) )
		f_new = itp_fl(w_new)/w_new
		iv_new = itp_iv(w_new)
		allflux[index] = f_new
		allivar[index] = iv_new


n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allflux_"+name+".txt")   , allflux)
n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allivar_"+name+".txt")   , allivar)
n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "masterwave_"+name+".txt"), masterwave)

