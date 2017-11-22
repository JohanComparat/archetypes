"""

comparat@ds52
cd  $HOME/SDSS/v5_10_7
scp -r u0936736@eboss.sdss.org:/uufs/chpc.utah.edu/common/home/sdss/ebosswork/eboss/spectro/redux/v5_10_7/spectra/full/9* .

#from importlib import reload
#import fitsio
#from SetCoverPy import setcover, mathutils
#import ebossspec
#import scipy.optimize as op
#from scipy.ndimage import gaussian_filter1d
#import cookb_signalsmooth

"""
import numpy as n
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
from SetCoverPy import setcover, mathutils

from scipy.stats import scoreatpercentile

import os
path_2_eboss21_cata = os.path.join(os.environ['OBS_REPO'], 'SDSS/ELG/eboss21/cats/eboss21.v5_10_7.latest.fits')

import astropy.io.fits as fits
cat = fits.open(path_2_eboss21_cata)[1].data

maskLambda = n.loadtxt(os.path.join(os.environ['GIT_ARCHETYPES'],'data',"dr12-sky-mask.txt"), unpack=True)

# parameters of the run
zmin = 0.7
zmax = 1.1
zQ_min = 3.
wlmin = 3600.
wlmax = 10300.
wlmin_rf = int(wlmin/(1+zmax))
wlmax_rf = int(wlmax/(1+zmin))
masterwave = n.arange(wlmin_rf, wlmax_rf, 0.5)

ok=(cat['Z_ZQ']>=zQ_min)&(cat['Z']>zmin)&(cat['Z']<=zmax)&(cat['CLASS']=="GALAXY")
#&(cat['PLATE']==9242)

galaxytype=cat['CLASS'][ok]
z=cat['Z'][ok]
zgood = n.ones_like(z)

nGal = len(ok.nonzero()[0])

name = "zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_zQmin_"+str(zQ_min)+"_Ngal_"+str(nGal)
print("run name:",name)
print("considering", nGal)

#allflux = n.loadtxt(os.path.join('..', 'archetypes_v0', "allflux.txt") )
#allivar = n.loadtxt(os.path.join('..', 'archetypes_v0', "allivar.txt") )
#masterwave = n.loadtxt(os.path.join('..', 'archetypes_v0', "masterwave.txt"))

allflux = n.loadtxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allflux_"+name+".txt") )
allivar = n.loadtxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "allivar_"+name+".txt") )
masterwave = n.loadtxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "masterwave_"+name+".txt"))

print("interpolates on the common grid")
index_wave_all = n.searchsorted(masterwave, [2400., 6200.])
tmpflux = allflux.T[index_wave_all[0]:index_wave_all[1],:]
tmpivar = allivar.T[index_wave_all[0]:index_wave_all[1],:]
tmpwave = masterwave[index_wave_all[0]:index_wave_all[1]]
tmploglam = n.log10(tmpwave)

#median_sn = n.zeros(nGal)
#for i in n.arange(nGal):
    #iuse = (n.where(tmpivar[:,i]>0))[0]
    #if iuse.size>0:
        #median_sn[i] = n.median(tmpflux[iuse,i]*n.sqrt(tmpivar[iuse,i]))

median_sn = n.median( tmpflux*tmpivar, axis=0)
 
iuse = n.where( median_sn>3.)[0]

# smooth the data
#n_pixels = len(tmpwave)
#newwave = n.median(tmpwave.reshape(n_pixels//4, 4), axis=1)
#newflux = n.sum(tmpflux.reshape(n_pixels//4, 4, tmpflux.shape[1]), axis=1)
#newivar = 1./n.sum(1./tmpivar.reshape(n_pixels//4, 4, tmpflux.shape[1]), axis=1)

newwave = tmpwave
newflux = tmpflux
newivar = tmpivar
tmpchi2 = n.zeros((iuse.size, iuse.size))
A = n.zeros((iuse.size, iuse.size))

print("creates the matrix")
tmp_yerr = 1./n.sqrt(newivar[:, iuse].T.reshape(iuse.size, newwave.size))
tmp_y = newflux[:,iuse].T
for i in n.arange(iuse.size):
    tmp_x = newflux[:, iuse[i]].T.reshape(1,newwave.size)
    tmp_xerr = 1./n.sqrt(newivar[:, iuse[i]].T.reshape(1,newwave.size))
    A_tmp, chi2_tmp = mathutils.quick_amplitude(tmp_x, tmp_y, tmp_xerr, tmp_yerr)
    A[i,:] = A_tmp
    tmpchi2[i,:] = chi2_tmp
    
chi2 = tmpchi2/(iuse.size-1) # reduced chi2

#pcs = n.array([1,10,25,50,75,90,99])
#pcs = n.array([16,17,18,19,20,21,22, 23, 24])
#scrs = scoreatpercentile(n.ravel(chi2), [16,17,18,19,20,21,22, 23, 24])

scr = scoreatpercentile(n.ravel(chi2), 20 )

chi2_min = scr # 0.05 # the minimum distance, the only free paramter
a_matrix = chi2<chi2_min # relationship matrix
cost = n.ones(iuse.size)
g = setcover.SetCover(a_matrix, cost)
# I'm using greedy just for demonstration
# g.greedy()
# SolveSCP() should be used to generate near-optimal solution
g.SolveSCP()
# These are the archetypes
iarchetype = n.nonzero(g.s)[0]
print(pcs[ii], scr, len(iarchetype), len(iarchetype)*1./nGal)

n_rep = n.sum(a_matrix[:, iarchetype], axis=0) # how many covered by the archetype?
isort = n.argsort(n_rep)[::-1]

tmpmedian = n.zeros((iarchetype.size, masterwave.size))

# These are the archetypal composites we want to use as the initial guess
for i in n.arange(iarchetype.size):
    itmp = a_matrix[:, iarchetype[i]] # These are the instances represented by the archetype
    for j in n.arange(masterwave.size):
        thisflux = allflux.T[j,iuse[itmp]]
        tmpmedian[i, j] = n.median(thisflux[thisflux!=0]) # Only use the objects that have this wavelength covered

#from _pickle import cPickle
#cPickle.dum

imax = n.max(isort)
print('imax', imax)
p.clf()
fig = p.figure(figsize=(10,imax*5))
fig.subplots_adjust(hspace=0)
for i in n.arange(0,imax,1):
	ax = fig.add_subplot(imax+1,1,i+1)
	ax.plot(masterwave[::3], tmpmedian[isort[i],:][::3])
	ax.set_xlim(2000, 8000)
	ax.set_ylim(-0.1, 3)
	ax.set_xticks([])
	ax.axvline(3727, color='k', ls='dashed')
	ax.axvline(5007, color='k', ls='dashed')
	print(i, n.count_nonzero(a_matrix[:,iarchetype[isort[i]]]))
	ax.grid()

#fig.set_xticks([2000, 3000, 4000, 5000])
fig.savefig(os.path.join(os.environ['OBS_REPO'], 'archetypes', "archetypes_"+name+".txt"))
#fig.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'sdss', 'elg', 'test.png'))
p.clf()

n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "archetypes_"+name+".txt")   , n.vstack((masterwave, tmpmedian)))
#n.savetxt(os.path.join('..', 'archetypes_v0', "archetypes.txt")   , n.vstack((masterwave, tmpmedian)))
# DATA = n.loadtxt(os.path.join('..', 'archetypes_v0', "archetypes.txt") )
