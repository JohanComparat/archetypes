import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import sys

# gets the eboss elg summary file
path_2_eboss21_cata = os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/spiders/target/2RXS_AllWISE_catalog_paper_2017May26_withSpectro_with2RXS_mask_DR14area_SPM_QSFIT_DR14Q.fits')
cat = fits.open(path_2_eboss21_cata)[1].data

# mask for the sky liens
maskLambda = n.loadtxt(os.path.join(os.environ['GIT_ARCHETYPES'],'data',"dr12-sky-mask.txt"), unpack=True)

get_path_to_eboss_spectrum = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', 'v5_10_0', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )
get_path_to_sdss_spectrum = lambda plate, mjd, fiberid : os.path.join(os.environ['HOME'], 'SDSS', '26', 'spectra', str(plate).zfill(4), "spec-"+str(plate).zfill(4)+"-"+str(mjd).zfill(5)+"-"+str(fiberid).zfill(4)+".fits" )

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
wlmax = 10300.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
masterwave = n.arange(wlmin_rf, wlmax_rf, 0.5)

ok=(cat['DR14_ZWARNING']==0)& (cat['DR14_Z']>zmin)& (cat['DR14_Z']<=zmax)

z=cat['DR14_Z'][ok]
zgood = n.ones_like(z)

nGal = len(ok.nonzero()[0])

name = "agn_zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_Ngal_"+str(nGal)
print("run name:",name)
print("considering", nGal)

allflux = n.zeros( ( nGal, len(masterwave)) )
allivar = n.zeros( ( nGal, len(masterwave)) )

print("gets the spectra to construct the matrix")
for index in range(nGal) :# = 0
	print(index)
	if cat['DR14_PLATE'][ok][index] < 3007:
	  path = get_path_to_sdss_spectrum(cat['DR14_PLATE'][ok][index], cat['DR14_MJD'][ok][index], cat['DR14_FIBERID'][ok][index])
	if cat['DR14_PLATE'][ok][index] > 3007:
	  path = get_path_to_eboss_spectrum(cat['DR14_PLATE'][ok][index], cat['DR14_MJD'][ok][index], cat['DR14_FIBERID'][ok][index])
	if os.path.isfile(path):
		print( path )
		redshift = cat['DR14_Z'][ok][index]
		wl, fl, iv = get_spec(path)
		#wl_rf = wl / (1+redshift)
		w_new = masterwave*(1+redshift)
		itp_fl = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,wl*fl,0.,0.)) )
		itp_iv = interp1d( n.hstack((100.,wl[0]-1.,wl,wl[-1]+1,20000.)), n.hstack((0.,0.,iv,0.,0.)) )
		f_new = itp_fl(w_new)/w_new
		iv_new = itp_iv(w_new)
		allflux[index] = f_new
		allivar[index] = iv_new
	else :
		print("not there", path )


n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allflux_"+name+".txt")   , allflux)
n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allivar_"+name+".txt")   , allivar)
n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "masterwave_"+name+".txt"), masterwave)
