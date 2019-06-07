import time
t0 = time.time()

import numpy as n
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as p
from scipy.interpolate import interp1d
from SetCoverPy import setcover, mathutils
import sys
from scipy.stats import scoreatpercentile
import os
import pickle
import SpectraStackingEBOSS as sse

version = 'v1'
sdss_dir = os.path.join(os.environ['HOME'], "SDSS")
stack_dir = os.path.join(sdss_dir, "stacks", "X_AGN")
archetype_dir = os.path.join(sdss_dir, "archetypes", version)

def stack_it(specList ):
 outfile = os.path.join(stack_dir, os.path.basename(specList)[:-4]+".stack")
 print(outfile)
 test_D = n.loadtxt(specList, unpack=True)
 print(len(test_D[0]))
 if len(test_D[0])>10:
  stack=sse.SpectraStackingEBOSS(specList, outfile )
  #stack.createStackMatrix()
  #stack.stackSpectra()
  return stack

file_list = n.array([
	os.path.join(stack_dir, "full_BLAGN_zmin_00_zmax_50.asc" ) ,
	os.path.join(stack_dir, "full_BLAZAR_zmin_00_zmax_50.asc") ,
	os.path.join(stack_dir, "full_BLLAC_zmin_00_zmax_50.asc" ) ,
	os.path.join(stack_dir, "full_GALAXY_zmin_00_zmax_50.asc") ,
	os.path.join(stack_dir, "full_NLAGN_zmin_00_zmax_50.asc" ) ,
	os.path.join(stack_dir, "full_NONE_zmin_00_zmax_50.asc"  ) ,
	os.path.join(stack_dir, "full_QSO_zmin_00_zmax_50.asc"   ) ,
	os.path.join(stack_dir, "full_STAR_zmin_00_zmax_50.asc"  ) ])


###############################################
###############################################
# NOW ARCHETYPE
###############################################
###############################################
class ObjIds:
	def __init__(self, imax, iuse):
		self.imax = imax
		self.iuse = iuse
		self.n_rep = n_rep

def make_archetype(stack, file_input, sn_min = 2.):
	print('starts archetype for',file_input, time.time()-t0)
	out_name = 'archetype_'+os.path.basename(stack[file_input].out_file)+'_snMin_'+str(sn_min)
	allflux    = n.loadtxt(stack[file_input].out_file+'.specMatrix.dat'      , unpack=True )      #, specMatrix)
	allivar    = n.loadtxt(stack[file_input].out_file+'.specMatrixErr.dat'   , unpack=True )**(-1)#, specMatrixErr)
	masterwave = stack[file_input].wave 
	tmploglam = n.log10(masterwave)
	# filter data
	median_sn = n.median( allflux*allivar, axis=0)
	iuse = n.where( median_sn>sn_min)[0]
	#
	tmpchi2 = n.zeros((iuse.size, iuse.size))
	A = n.zeros((iuse.size, iuse.size))
	print("creates the matrix", time.time()-t0, 's', allflux.shape)
	tmp_yerr = 1./n.sqrt(allivar[:, iuse].T.reshape(iuse.size, masterwave.size))
	tmp_y = allflux[:,iuse].T
	for i in n.arange(iuse.size):
		#print(i)
		tmp_x = allflux[:, iuse[i]].T.reshape(1,masterwave.size)
		tmp_xerr = 1./n.sqrt(allivar[:, iuse[i]].T.reshape(1,masterwave.size))
		#print(tmp_x, tmp_y, tmp_xerr, tmp_yerr)
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

	# sort archeypes against the number of spectra they represent
	n_rep = n.sum(a_matrix[:, iarchetype], axis=0) # how many covered by the archetype?
	isort = n.argsort(n_rep)[::-1]
	print(n_rep)
	archetype_median = n.zeros((iarchetype.size, masterwave.size))

	# These are the archetypal composites we want to use as the initial guess
	for i in n.arange(iarchetype.size):
		itmp = a_matrix[:, iarchetype[i]] # These are the instances represented by the archetype
		for j in n.arange(masterwave.size):
			thisflux = allflux[j,iuse[itmp]]
			archetype_median[i, j] = n.median(thisflux[thisflux!=0]) # Only use the objects that have this wavelength covered

	#from _pickle import cPickle
	#cPickle.dum

	imax = n.max(isort)
	print('imax', imax)
	p.clf()
	fig = p.figure(figsize=(10,imax*5))
	fig.subplots_adjust(hspace=0)
	for i in n.arange(0,imax,1):
		ax = fig.add_subplot(imax+1,1,i+1)
		ax.plot(masterwave[::3], archetype_median[isort[i],:][::3] , label = 'nRep=' + str(n_rep[isort[i]]) )
		ax.set_xlim(masterwave[10], masterwave[-10])
		#ax.set_ylim(-0.1, 3)
		ax.set_xticks([])
		ax.axvline(1215, color='b', ls='dashed', label='1215 Lya')
		ax.axvline(1546, color='c', ls='dashed', label='1546 CIV')
		ax.axvline(2800, color='m', ls='dashed', label='2800 MgII')
		ax.axvline(3727, color='g', ls='dashed', label='3727 [OII]')
		ax.axvline(5007, color='r', ls='dashed', label='5007 [OIII]')
		ax.axvline(6565, color='k', ls='dashed', label='6565 Ha')
		print(i, n.count_nonzero(a_matrix[:,iarchetype[isort[i]]]))
		ax.grid()
		ax.legend(frameon=False, loc=0)


	fig.savefig( os.path.join(archetype_dir, "figure_archetypes_"+out_name +".png") )
	p.clf()

	n.savetxt(os.path.join(archetype_dir, "archetypes_"+out_name+".txt")   , n.vstack((masterwave, archetype_median)))
	
	f = open(os.path.join(archetype_dir, "index_archetypes_"+out_name+".pkl"), 'wb')
	obj = ObjIds(imax, iuse, n_rep)
	pickle.dump(obj, f)
	f.close()


stack = {}
#for file_input in file_list:
file_input = file_list[2]
stack[file_input] = stack_it(file_input)
make_archetype(stack, file_input, sn_min = 2.)
#n.savetxt(os.path.join('..', 'archetypes_v0', "archetypes.txt")   , n.vstack((masterwave, archetype_median)))
# DATA = n.loadtxt(os.path.join('..', 'archetypes_v0', "archetypes.txt") )
