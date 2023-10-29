#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 09:53:18 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import subprocess, os, shutil
band ='ir'
img = '/scratch/bell/dutta26/backup/wted_coadds/ir_coadd_wt.fits'
sextractorFile = '/scratch/bell/dutta26/wiyn_sim/coadd_sim.cat'
starGal_out_file = '/scratch/bell/dutta26/backup/source_list.pk1'
coadd_detect_file = '/scratch/bell/dutta26/backup/coaddSc_' + str(band)+'.npy'
ebModeFile = '/scratch/bell/dutta26/backup/master_arr_coaddMC.npy'
plotLoc = '/scratch/bell/dutta26/backup/plot/'
bandLoc = '/scratch/bell/dutta26/backup/'
ir_band_data = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'


psfMin = 3.3
psfMax = 3.8
fluxMin = 5e2
fluxMax = 1e5

psfMin_2 = 2.5
psfMax_2 = 4.5
fluxMin_2 = 1e2
fluxMax_2 = 1e6

# =============================================================================
# #Run sextractor
# sexLoc = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/'
# sexCommand = 'sex -c /scratch/bell/dutta26/wiyn_sim/default.sex '+str(med_img)+' -DETECT_MINAREA 5 -DETECT_THRESH 0.8 -ANALYSIS_THRESH 0.8 -WEIGHT_TYPE BACKGROUND -CATALOG_NAME '+sextractorFile
# process = subprocess.Popen(sexCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
# output, error = process.communicate()
# 
# =============================================================================
# =============================================================================
# import seperate_star_gal
# seperate_star_gal.seperate(sextractorFile, img, starGal_out_file, plotLoc, psfMin, psfMax, fluxMin, fluxMax, psfMin_2, psfMax_2, fluxMin_2, fluxMax_2)
# =============================================================================

# =============================================================================
# import coadd_detect
# coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)
# =============================================================================

import makeEandBmodes
makeEandBmodes.EandB( ir_band_data,  '/home/dutta26/codes/wiyn_wl_sim/coaddSc_r.npy', '/home/dutta26/codes/wiyn_wl_sim/coaddSc_i.npy', None, ebModeFile, img, '/scratch/bell/dutta26/wiyn_sim/r_withoutMC.npy', '/scratch/bell/dutta26/wiyn_sim/i_withoutMC.npy')





# =============================================================================
# band ='u'
# img = '/scratch/bell/dutta26/backup/wted_coadds/'+band+'_coadd_wt.fits'
# sextractorFile = '/scratch/bell/dutta26/wiyn_sim/coadd_sim.cat'
# starGal_out_file = '/scratch/bell/dutta26/backup/source_list.pk1'
# coadd_detect_file = '/scratch/bell/dutta26/backup/coaddSc_' + str(band)+'.npy'
# plotLoc = '/scratch/bell/dutta26/backup/plot/'
# bandLoc = '/scratch/bell/dutta26/backup/'
# ir_band_data = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'
# 
# import coadd_detect
# coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)
# 
# 
# =============================================================================
# =============================================================================
# band ='g'
# img = '/scratch/bell/dutta26/backup/wted_coadds/'+band+'_coadd_wt.fits'
# sextractorFile = '/scratch/bell/dutta26/wiyn_sim/coadd_sim.cat'
# starGal_out_file = '/scratch/bell/dutta26/backup/source_list.pk1'
# coadd_detect_file = '/scratch/bell/dutta26/backup/coaddSc_' + str(band)+'.npy'
# plotLoc = '/scratch/bell/dutta26/backup/plot/'
# bandLoc = '/scratch/bell/dutta26/backup/'
# ir_band_data = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'
# 
# 
# import coadd_detect
# coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)
# 
# =============================================================================




