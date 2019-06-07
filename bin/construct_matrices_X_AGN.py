"""
python3.4 construct_matrices_sdss.py -0.1 0.05 X_AGN
"""
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle
import spectres
import SpectraStackingEBOSS as sse

version = 'v1'
sdss_dir = = os.path.join(os.environ['HOME'], "SDSS")
archetype_dir = os.path.join(sdss_dir, "archetypes", version)

logwlmin = 2.9. #1000.
logwlmax = 4.0. #10000.
d_lambda = 0.0001 
masterwave = 10**n.arange(logwlmin, logwlmax, d_lambda) 


file_list = n.array([
	"/home/comparat/SDSS/stacks/X_AGN/full_BLAGN_zmin_00_zmax_50.asc"  ,
	#"/home/comparat/SDSS/stacks/X_AGN/full_BLANK_zmin_00_zmax_50.asc"  ,
	"/home/comparat/SDSS/stacks/X_AGN/full_BLAZAR_zmin_00_zmax_50.asc" ,
	"/home/comparat/SDSS/stacks/X_AGN/full_BLLAC_zmin_00_zmax_50.asc"  ,
	"/home/comparat/SDSS/stacks/X_AGN/full_GALAXY_zmin_00_zmax_50.asc" ,
	"/home/comparat/SDSS/stacks/X_AGN/full_NLAGN_zmin_00_zmax_50.asc"  ,
	"/home/comparat/SDSS/stacks/X_AGN/full_NONE_zmin_00_zmax_50.asc"   ,
	"/home/comparat/SDSS/stacks/X_AGN/full_QSO_zmin_00_zmax_50.asc"    ,
	"/home/comparat/SDSS/stacks/X_AGN/full_STAR_zmin_00_zmax_50.asc"   ])


def stack_it(specList ):
	outfile = join(spec_dir, os.path.basename(specList)[:-4]+".stack")
	print(outfile)
	test_D = n.loadtxt(specList, unpack=True)
	print(len(test_D[0]))
	if len(test_D[0])>10:
		stack=sse.SpectraStackingEBOSS(specList, outfile )
		stack.createStackMatrix()

		#stack.stackSpectra()

for file_input in file_list[4:]:
	stack_it(file_input)


def create_matrix(path_2_list):
path_2_list = all_lists[0]


nGal = len(ok.nonzero()[0])
# 1. for each element of the list get the spectrum
# 2. in the rest frame force it to the masterwave

print("run name:",name)
print("considering", nGal)
SSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
allflux = n.zeros( ( nGal, len(masterwave)) )v5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
allivar = n.zeros( ( nGal, len(masterwave)) )

print("gets the spectra to construct the matrix")
for index in range(nGal) :# = 0
	path_1 = get_path_to_spectrum_v5_10_10(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	path_2 = get_path_to_spectrum_26(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
	path_3 = get_path_to_spectrum_v5_10_0(cat['PLATE'][ok][index], cat['MJD'][ok][index], cat['FIBERID'][ok][index])
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
		max_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))
out_dir = os.path.join( os.environ['OBS_REPO'], 'archetypes', version, gal_type )
if os.path.isdir(out_dir)==False:], cat['FIBERID_BEST'][ok_i], z_best[ok_i] ])
  os.system('mkdir -p '+out_dir)
  
n.savetxt(os.path.join(out_dir, "allflux_"+name+".txt")   , allflux)
n.savetxt(os.path.join(out_dir, "allivar_"+name+".txt")   , allivar)
n.savetxt(os.path.join(out_dir, "masterwave_"+name+".txt"), masterwave)                                                                                                            


class ObjIds:                                                                                                                  
    def __init__(self, plate, mjd, fiberid): QSO_BAL'], dtype='<U16')
        self.plate = plate
        self.mjd = mjd                                                                                                                   
        self.fiberid = fiberid
        
f = open(os.path.join(out_dir, "specObj_"+name+".pkl"), 'wb')
obj = ObjIds(cat['PLATE'][ok], cat['MJD'][ok], cat['FIBERID'][ok])
pickle.dump(obj, f)
f.close()


	wavelength_spectrum = np.hstack(( 
		wave[0]-10, 
		wave[0]-5,
		np.min(inter_lambda)-10,
		np.min(inter_lambda)-5,
		inter_lambda,
		np.max(inter_lambda)+5,
		np.max(inter_lambda)+10,
		wave[-1]+5, 
		wave[-1]+10
		))
	#
	flux_spectrum = np.hstack(( 
		inter_flux[0],inter_flux[0],inter_flux[0],inter_flux[0],
		inter_flux,
		inter_flux[-1],inter_flux[-1],inter_flux[-1],inter_flux[-1]
		))
	#
	print(flux_spectrum.shape)
	flux_error_spectrum = flux_spectrum/100.
	#
	final_spectrum, final_spectrum_err = sp.spectres(
		wave, 
		wavelength_spectrum, 
		flux_spectrum, 
		flux_error_spectrum )




