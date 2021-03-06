"""
exploitation script: plots archetypes looks for the interesting ones in the data.
"""
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

version = 'v6'

# parameters of the run
zmin = float(sys.argv[1])
zmax = float(sys.argv[2])
gal_type = sys.argv[3]
# ELG, LRG, X_AGN, QSO

Nspec_max = 3000

name = "sdss_"+gal_type+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))

print(gal_type)
print(name)
sn_min = float(sys.argv[4])

zmin = float(zmin)
zmax = float(zmax)

out_dir = os.path.join( os.environ['OBS_REPO'], 'archetypes', version, gal_type )


allflux = n.loadtxt(os.path.join(out_dir, "allflux_"+name+".txt") )
allivar = n.loadtxt(os.path.join(out_dir, "allivar_"+name+".txt") )
masterwave = n.loadtxt(os.path.join(out_dir, "masterwave_"+name+".txt"))

#  , n.vstack((masterwave, tmpmedian))
DATA = n.loadtxt(os.path.join(out_dir, "archetypes_"+name+"_snMin"+str(sn_min)+".txt")  )


imax = len(DATA)
fig = p.figure(figsize=(10,imax*5))
fig.subplots_adjust(hspace=0)
for i in n.arange(1,imax,1):
	ax = fig.add_subplot(imax-1,1,i)
	ax.plot(DATA[0], DATA[i]/n.median(DATA[i]) )
	ax.set_ylim(0.5, 10)
	ax.set_yscale('log')
	ax.set_xticks([])
	ax.axvline(1215, color='b', ls='dashed', label='1215 Lya')
	ax.axvline(1546, color='c', ls='dashed', label='1546 CIV')
	ax.axvline(2800, color='m', ls='dashed', label='2800 MgII')
	ax.axvline(3727, color='g', ls='dashed', label='3727 [OII]')
	ax.axvline(5007, color='r', ls='dashed', label='5007 [OIII]')
	ax.axvline(6565, color='k', ls='dashed', label='6565 Ha')
	#print(i, n.count_nonzero(a_matrix[:,iarchetype[isort[i]]]))
	ax.grid()
	ax.set_xlim((DATA[0].min(), DATA[0].max()))

ax.legend(frameon=False)
ax.set_xlabel('wavelength')
fig.savefig(os.path.join(out_dir, "archetypes_"+name+"_snMin"+str(sn_min)+".png"))
p.clf()

sys.exit()

#f = open(os.path.join(out_dir, "specObj_"+name+".pkl"), 'rb')
#obj_in = pickle.load(f)
#f.close()


#f = open(os.path.join(out_dir, "indexArchetypes_"+name+"_snMin"+str(sn_min)+".pkl"), 'r')
#rep = pickle.load( f)
#f.close()

#sys.exit()

#import time
#t0 = time.time()

#import numpy as n
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as p

#import pickle


print('imax', rep.imax)
p.clf()
fig = p.figure(figsize=(10,rep.imax*5))
fig.subplots_adjust(hspace=0)
for i in n.arange(0,rep.imax,1):
	ax = fig.add_subplot(rep.imax+1,1,i+1)
	ax.plot(masterwave[::3], tmpmedian[isort[i],:][::3] , label = 'nSpec=' + str(n.count_nonzero(a_matrix[:,iarchetype[isort[i]]])) )
	ax.set_xlim(2000, 8000)
	ax.set_ylim(-0.1, 3)
	ax.set_xticks([])
	ax.axvline(1215, color='b', ls='dashed', label='1215 Lya')
	ax.axvline(1546, color='c', ls='dashed', label='1546 CIV')
	ax.axvline(2800, color='m', ls='dashed', label='2800 MgII')
	ax.axvline(3727, color='g', ls='dashed', label='3727 [OII]')
	ax.axvline(5007, color='r', ls='dashed', label='5007 [OIII]')
	ax.axvline(6565, color='k', ls='dashed', label='6565 Ha')
	print(i, n.count_nonzero(a_matrix[:,iarchetype[isort[i]]]))
	ax.grid()

#fig.set_xticks([2000, 3000, 4000, 5000])
#fig.savefig(os.path.join(os.environ['OBS_REPO'], 'archetypes', "archetypes_"+name+"_snMin"+str(sn_min)+".png"))
#fig.savefig(os.path.join(os.environ['HOME'], 'wwwDir', 'sdss', 'elg', 'test.png'))
p.clf()

n.savetxt(os.path.join(os.environ['OBS_REPO'], 'archetypes', "archetypes_"+name+"_snMin"+str(sn_min)+".txt")   , n.vstack((masterwave, tmpmedian)))
#n.savetxt(os.path.join('..', 'archetypes_v0', "archetypes.txt")   , n.vstack((masterwave, tmpmedian)))
# DATA = n.loadtxt(os.path.join('..', 'archetypes_v0', "archetypes.txt") )
