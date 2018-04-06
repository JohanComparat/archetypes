"""
Routine to create the lists of spectra to be stacked

Then stacks them using the class :
/home/comparat/software/linux/pySU/galaxy/python/SpectraStackingEBOSS.py


"""
import os
import sys
import SpectraStackingEBOSS as sse

gal_type= sys.argv[1]
zmin = float(sys.argv[2])
zmax = float(sys.argv[3])
#gal_type='qso_BL'
#zmin = 0.
#zmax = 0.5
Nspec_max = 10000

version = 'v1'

name = "sdss_"+gal_type+"_zmin_"+str(int(10*zmin)).zfill(2)+"_zmax_"+str(int(10*zmax)).zfill(2)+"_Nlt_"+str(int(Nspec_max))


out_dir = os.path.join( os.environ['OBS_REPO'], 'spectrastacks', version, gal_type )
list_2_stack = os.path.join(out_dir,  "stacklist_"+name+".ascii")
output_file_name = os.path.join(out_dir,  "spec_"+name+".fits")


st=sse.SpectraStackingEBOSS(list_2_stack, output_file_name)
st.stackSpectra()

