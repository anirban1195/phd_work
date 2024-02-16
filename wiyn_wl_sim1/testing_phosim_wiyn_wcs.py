#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 16:57:12 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np
from astropy.io import fits
import helper_phosim
import subprocess,sys
phosimLoc = '/home/dutta26/Downloads/phosim_release/'
import gzip,shutil


phosim_command = './phosim /scratch/bell/dutta26/wiyn_sim/test_star.txt -i wiyn_odi -c examples/nobackground -o /scratch/bell/dutta26/wiyn_sim/test/ -w /scratch/bell/dutta26/wiyn_sim/work/ -p 5 -t 5 '
process = subprocess.Popen(phosim_command.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
output, error = process.communicate()

phosim_command = './phosim /scratch/bell/dutta26/wiyn_sim/test_star.txt -i generic_2m -c examples/nobackground -o /scratch/bell/dutta26/wiyn_sim/test/ -w /scratch/bell/dutta26/wiyn_sim/work/ -p 5 -t 5 '
process = subprocess.Popen(phosim_command.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
output, error = process.communicate()

with gzip.open('/scratch/bell/dutta26/wiyn_sim/test/generic_2m_e_1098_f2_chip_E000.fits.gz', 'rb') as f_in:
    with open('/scratch/bell/dutta26/wiyn_sim/test/generic_2m_e_1098_f2_chip_E000.fits', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        
with gzip.open('/scratch/bell/dutta26/wiyn_sim/test/wiyn_odi_e_1098_f2_OTA44_E000.fits.gz', 'rb') as f_in:
    with open('/scratch/bell/dutta26/wiyn_sim/test/wiyn_odi_e_1098_f2_OTA44_E000.fits', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)