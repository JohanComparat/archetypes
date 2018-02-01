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
z_all = n.hstack(( 0.001, n.arange(0.1, 3., 0.2) ))


command_plot = lambda zmin, zmax, gal_type : "python3.4 plot_archetypes.py "+str(zmin)+" "+str(zmax)+" "+gal_type+" 3."

for zmin, zmax in zip(z_all[:-1], z_all[1:]):
	cm_plot = command_plot(zmin, zmax, gal_tp)
	print(cm_plot)
	os.system(cm_plot)
