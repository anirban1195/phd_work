#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 23:18:29 2023

@author: dutta26
"""

import numpy as np
from astropy.io import fits
import helper_phosim
import subprocess,sys

phosimLoc = '/home/dutta26/Downloads/phosim_release/'

filt = str(sys.argv[1])
idNo = 1000
noImages = 30

for j in range(noImages):
    catalog = helper_phosim.make_phosim_imputCatalog(str(idNo), filt)
    catalog_star = helper_phosim.make_phosim_star_imputCatalog(str(idNo), filt)
    
    #For normal images
    outLoc = '/scratch/bell/dutta26/wiyn_sim/'+filt+'/'+str(idNo)+"/"
    helper_phosim.check_n_remove_existing(outLoc)
    phosim_command = './phosim '+catalog+' -i wiyn_odi -c examples/quickbackground -o '+outLoc+' -w /scratch/bell/dutta26/wiyn_sim/work1/ -s "OTA32|OTA33|OTA34|OTA35|OTA43|OTA44|OTA45|OTA53|OTA54|OTA55|OTA23" '
    #phosim_command = './phosim '+catalog+' -i generic_4m -c examples/quickbackground -o '+outLoc+' -w /scratch/bell/dutta26/wiyn_sim/work/'
    process = subprocess.Popen(phosim_command.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
    output, error = process.communicate()
    helper_phosim.clean(outLoc, catalog)
    del catalog, outLoc
    
    
    
    #For star images
    outLoc_star = '/scratch/bell/dutta26/wiyn_sim/'+filt+'/'+str(idNo)+"_star/"
    helper_phosim.check_n_remove_existing(outLoc_star)
    phosim_command = './phosim '+catalog_star+' -i wiyn_odi -c examples/quickbackground -o '+outLoc_star+' -w /scratch/bell/dutta26/wiyn_sim/work1/'
    #phosim_command = './phosim '+catalog_star+' -i generic_2m -c examples/nobackground -o '+outLoc_star+' -w /scratch/bell/dutta26/wiyn_sim/work/'
    process = subprocess.Popen(phosim_command.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
    output, error = process.communicate()
    helper_phosim.clean(outLoc_star, catalog_star)
    del outLoc_star, catalog_star
    
    
    idNo += 1