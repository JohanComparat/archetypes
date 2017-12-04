#!/bin/bash

python3.4 construct_matrices_sdss_lrg_spectra.py 0.0 0.1
python3.4 construct_matrices_sdss_lrg_spectra.py 0.1 0.2
python3.4 construct_matrices_sdss_lrg_spectra.py 0.2 0.3
python3.4 construct_matrices_sdss_lrg_spectra.py 0.3 0.4
python3.4 construct_matrices_sdss_lrg_spectra.py 0.4 0.5
python3.4 construct_matrices_sdss_lrg_spectra.py 0.5 0.6
python3.4 construct_matrices_sdss_lrg_spectra.py 0.6 0.7
python3.4 construct_matrices_sdss_lrg_spectra.py 0.7 0.8
python3.4 construct_matrices_sdss_lrg_spectra.py 0.8 0.9
python3.4 construct_matrices_sdss_lrg_spectra.py 0.9 1.0
python3.4 construct_matrices_sdss_lrg_spectra.py 1.0 1.2

python3.4 construct_archetypes.py ebossLRG_zmin_0.0_zmax_0.1_Ngal_5041 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.1_zmax_0.2_Ngal_4964 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.2_zmax_0.3_Ngal_4953 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.3_zmax_0.4_Ngal_5022 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.4_zmax_0.5_Ngal_5005 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.5_zmax_0.6_Ngal_5041 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.6_zmax_0.7_Ngal_4906 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.7_zmax_0.8_Ngal_5045 2.
python3.4 construct_archetypes.py ebossLRG_zmin_0.8_zmax_0.9_Ngal_4980 2.




