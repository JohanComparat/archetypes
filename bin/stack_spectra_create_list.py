"""
Routine to create the lists of spectra to be stacked

Then stacks them using the class :
/home/comparat/software/linux/pySU/galaxy/python/SpectraStackingEBOSS.py


"""
import numpy as n
import astropy.io.fits as fits
import sys
import os

gal_type= sys.argv[1]
zmin = float(sys.argv[2])
zmax = float(sys.argv[3])
#gal_type='qso_BL'
#zmin = 0.
#zmax = 0.5
Nspec_max = 10000

version = 'v1'

name = "sdss_"+gal_type+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))

output_file = "stacklist_"+name+".ascii"

out_dir = os.path.join( os.environ['OBS_REPO'], 'spectrastacks', version, gal_type )
if os.path.isdir(out_dir)==False:
  os.system('mkdir -p '+out_dir)

path_2_outfile = os.path.join(out_dir, output_file)

if gal_type=='qso_BL':
  cat = fits.open('/home/comparat/SDSS/catalogs/DR14Q_v4_4.fits')[1].data
  ok_i=(cat['Z']>zmin)&(cat['Z']<=zmax)
  print(len(cat['Z'][ok_i]))
  #print(n.histogram(cat['Z'][ok],bins=n.arange(0,3.1,0.5)))


if gal_type == "LRG":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/dr14/specObj-dr14.fits'))[1].data
	ok_i=(cat['Z']>zmin)&(cat['Z']<=zmax)&((cat['SOURCETYPE']=='CMASS')|(cat['SOURCETYPE']=='LRG')|(cat['SOURCETYPE']=='FAINT_HIZ_LRG')|(cat['SOURCETYPE']=='HIZ_LRG')|(cat['SOURCETYPE']=='LRG_ROUND3'))
	print(len(cat['Z'][ok_i]))
	print(n.histogram(cat['Z'][ok_i],bins=n.arange(0,3.1,0.5)))

if gal_type == "ELG":
	cat = fits.open(os.path.join(os.environ['OBS_REPO'], 'SDSS/ELG/eboss21/cats/eboss21.v5_10_7.latest.fits'))[1].data
	zQ_min = 2.
	ok_i=(cat['Z_ZQ']>=zQ_min)&(cat['Z']>zmin)&(cat['Z']<=zmax)&(cat['SN_MEDIAN_ALL']>0.)
	print(len(cat['Z'][ok_i]))
	print(n.histogram(cat['Z'][ok_i],bins=n.arange(0,3.1,0.5)))


nGal_all = len(ok_i.nonzero()[0])
print("total", nGal_all)
if nGal_all > Nspec_max:
    rds = n.random.rand(len(cat))
    ok = (ok_i)&(rds < Nspec_max / nGal_all)
else :
    ok = ok_i
    
nGal = len(ok.nonzero()[0])
print("to the list", nGal)

n.savetxt(path_2_outfile, n.transpose([cat['PLATE'][ok], cat['MJD'][ok], cat['FIBERID'][ok], cat['Z'][ok]]) )


