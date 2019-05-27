import os
import sys
import numpy as n

gal_tp = sys.argv[1]

#qso_BL
#qso_t2
#X_no_BL
#GAL_agn_all
#XAGN
#XMMSL
#QSO
#ELG
#LRG

#z_all = n.hstack((-0.001, 0.001, n.arange(0.1, 3., 0.2) ))
#z_all = n.hstack(( 0.001, n.arange(0.1, 3., 0.2) ))
z_all=n.array([3., 6.])

command = lambda zmin, zmax, gal_type : "python3.4 construct_matrices_sdss.py "+str(zmin)+" "+str(zmax)+" "+gal_type

for zmin, zmax in zip(z_all[:-1], z_all[1:]):
	cm = command(zmin, zmax, gal_tp)
	print(cm)
	os.system(cm)
