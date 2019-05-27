#! /usr/bin/env python

"""
This script produces the stacks for emission line luminosity limited samples.

nohup python2 stack_spectra_ELG_LineLF.py > stack_spectra_ELG_LineLF.log &

"""
import sys
import os 
from os.path import join
import glob
import numpy as n
import SpectraStackingEBOSS as sse

spec_dir = join(os.environ['HOME'],"SDSS/stacks/X_AGN")

file_list = n.array(glob.glob(os.path.join(spec_dir, 'full_*.ascii')))

def stack_it(specList ):
	outfile = join(spec_dir, os.path.basename(specList)[:-6]+".stack")
	stack=sse.SpectraStackingEBOSS(specList, outfile   )
	stack.createStackMatrix_Weighted()
	print(outfile)
	stack.stackSpectra()

for file_input in file_list:
	stack_it(file_input)


