"""
python3.4 construct_matrices_sdss.py -0.1 0.05 X_AGN
"""
import sys
import numpy as n
from scipy.interpolate import interp1d
import os
import astropy.io.fits as fits
import pickle

version = 'v3'

# parameters of the run
zmins = n.arange(0,3.1,0.5)[:-1]
zmaxs = n.arange(0,3.1,0.5)[1:]
#gal_type = sys.argv[3]
# ELG, LRG, X_AGN, QSO
#print('sys.argv', sys.argv)


SN_threshold = 0.5
wlmin = 3700.
wlmax = 9000.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
d_lambda = 1. 
masterwave = n.arange(wlmin_rf, wlmax_rf, d_lambda) #0.5)

path_2_cats = os.path.join(os.environ['HOME'], 'data', 'spiders', 'agn')
path_2_stack_lists = os.path.join(os.environ['HOME'], 'SDSS/stacks/X_AGN')


path_2_XMMSL2 = os.path.join(path_2_cats, 'SPIDERS_XMMSL2_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
path_2_2RXS   = os.path.join(path_2_cats, 'SPIDERS_2RXS_Xray_NWAY_ALLWISE_SDSSv5b_SpecDR14_SpecDR16_with_VI_1rowperXray_inDR16wSEQUELS_COMPLETE.fits')
path_2_XXL    = os.path.join(path_2_cats, 'Menzel16_spAll_XXL_plate.fits')
path_2_S82X   = os.path.join(path_2_cats, 'LaMassa19_spAll_S82X_plate.fits')

cat_XMMSL2 = fits.open(path_2_XMMSL2)[1].data
cat_2RXS   = fits.open(path_2_2RXS  )[1].data
cat_XXL    = fits.open(path_2_XXL   )[1].data
cat_S82X   = fits.open(path_2_S82X  )[1].data

# join the 4 catalogs into a single list for each category

# XMMSL2 lists
survey = 'XMMSL2'
cat = cat_XMMSL2
z_best = cat['Z_BEST']
z_conf = cat['CONF_BEST']
print(survey)
for zmin, zmax in zip(zmins, zmaxs):
	print('=================================================')
	print(zmin, '<z<', zmax)
	selection = (z_best > zmin) & (z_best < zmax) & (z_conf>=3)
	all_categories = n.unique(cat['CLASS_BEST'][selection])              
	for category in all_categories:
		ok_i = (cat['CLASS_BEST']==category) & (selection)
		N_obj = len(ok_i.nonzero()[0])
		print(category, N_obj)
		if N_obj>0:
			name = survey+"_"+category+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))
			print(name)
			DATA = n.transpose([ cat['PLATE_BEST'][ok_i], cat['MJD_BEST'][ok_i], cat['FIBERID_BEST'][ok_i], z_best[ok_i] ])
			n.savetxt(os.path.join(path_2_stack_lists, name), DATA)

sys.exit()

n.unique(cat_2RXS['CLASS_BEST'])                                                                                                                                                          


n.unique(cat_XXL['TDVI_CLASS_PERSON'])                                                                                                                                                          
# chararray(['', '     QSO', '    NONE', '    STAR', '  BLAZAR', '  GALAXY', ' QSO_BAL'], dtype='<U16')

n.unique(cat_S82X['TDVI_CLASS_PERSON'])                                                                                                                                                          
# chararray(['', '     QSO', '    NONE', '  GALAXY'], dtype='<U16')

cat = cat_XMMSL2
print(cat['PLATE_BEST'][:10], cat['MJD_BEST'][:10], cat['FIBERID_BEST'][:10])

cat = cat_2RXS
print(cat['PLATE_BEST'][:10], cat['MJD_BEST'][:10], cat['FIBERID_BEST'][:10])

cat = cat_XXL
print(cat['DR16_PLATE'][:10], cat['DR16_MJD'][:10], cat['DR16_FIBERID'][:10])

cat = cat_S82X
print(cat['DR16_PLATE'][:10], cat['DR16_MJD'][:10], cat['DR16_FIBERID'][:10])
