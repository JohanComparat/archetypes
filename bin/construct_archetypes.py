import time
t0 = time.time()

import numpy as n
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
from SetCoverPy import setcover, mathutils
import sys
from scipy.stats import scoreatpercentile
import os
import pickle

version = 'v3'

# parameters of the run
zmin = sys.argv[1]
zmax = sys.argv[2]
gal_type = sys.argv[3]
# ELG, LRG, X_AGN, QSO
Nspec_max = 2000.
name = "sdss_"+gal_type+"_zmin_"+str(zmin)+"_zmax_"+str(zmax)+"_Nlt_"+str(Nspec_max)

zmin = float(zmin)
zmax = float(zmax)

print(gal_type)
print(name)
sn_min = float(sys.argv[4])



out_dir = os.path.join( os.environ['OBS_REPO'], 'archetypes', version, gal_type )


allflux = n.loadtxt(os.path.join(out_dir, "allflux_"+name+".txt") )
allivar = n.loadtxt(os.path.join(out_dir, "allivar_"+name+".txt") )
masterwave = n.loadtxt(os.path.join(out_dir, "masterwave_"+name+".txt"))

print("interpolates on the common grid", time.time()-t0, 's')
index_wave_all = n.searchsorted(masterwave, [masterwave[0]+100,masterwave[-1]-100])
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
 
iuse = n.where( median_sn>sn_min)[0]

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

print("creates the matrix", time.time()-t0, 's')
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
print( scr, len(iarchetype), time.time()-t0, 's')

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
	ax.plot(masterwave[::3], tmpmedian[isort[i],:][::3] , label = 'nSpec=' + str(n.count_nonzero(a_matrix[:,iarchetype[isort[i]]])) )
	ax.set_xlim(masterwave[10], masterwave[-10])
	#ax.set_ylim(-0.1, 3)
	ax.set_xticks([])
	ax.axvline(1215, color='b', ls='dashed', label='1215')
	ax.axvline(2800, color='m', ls='dashed', label='2800')
	ax.axvline(3727, color='g', ls='dashed', label='3727')
	ax.axvline(5007, color='r', ls='dashed', label='5007')
	print(i, n.count_nonzero(a_matrix[:,iarchetype[isort[i]]]))
	ax.grid()
	ax.legend(frameon=False, loc=0)

#fig.set_xticks([2000, 3000, 4000, 5000])
fig.savefig(os.path.join(out_dir, "figureArchetypes_"+name+"_snMin"+str(sn_min)+".png"))
#fig.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'sdss', 'elg', 'test.png'))
p.clf()

n.savetxt(os.path.join(out_dir, "archetypes_"+name+"_snMin"+str(sn_min)+".txt")   , n.vstack((masterwave, tmpmedian)))

class ObjIds:
    def __init__(self, imax, iuse):
        self.imax = imax
        self.iuse = iuse
        
f = open(os.path.join(out_dir, "indexArchetypes_"+name+"_snMin"+str(sn_min)+".pkl"), 'wb')
obj = ObjIds(imax, iuse)
pickle.dump(obj, f)
f.close()

#n.savetxt(os.path.join('..', 'archetypes_v0', "archetypes.txt")   , n.vstack((masterwave, tmpmedian)))
# DATA = n.loadtxt(os.path.join('..', 'archetypes_v0', "archetypes.txt") )
