#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:19:32 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import subprocess, os, shutil
band ='ir'
img = '/scratch/bell/dutta26/wiyn_sim1/wted_coadds/ir_coadd_wt.fits'
med_img = '/scratch/bell/dutta26/wiyn_sim1/median_coadds/ir_coadd_med_crop.fits'
sextractorFile = '/scratch/bell/dutta26/wiyn_sim1/coadd_sim.cat'
starGal_out_file = '/home/dutta26/codes/wiyn_wl_sim1/source_list.pk1'
coadd_detect_file = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_' + str(band)+'.npy'
ebModeFile = '/scratch/bell/dutta26/wiyn_sim1/master_arr_coaddMC_sfMC.npy'
plotLoc = '/scratch/bell/dutta26/wiyn_sim1/plot/phosim/'
bandLoc = '/scratch/bell/dutta26/wiyn_sim1/'
ir_band_data = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy'


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
# sexCommand = 'sex -c /scratch/bell/dutta26/wiyn_sim1/default.sex '+str(med_img)+' -DETECT_MINAREA 5 -DETECT_THRESH 0.8 -ANALYSIS_THRESH 0.8 -WEIGHT_TYPE BACKGROUND -CATALOG_NAME '+sextractorFile
# process = subprocess.Popen(sexCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
# output, error = process.communicate()
# =============================================================================

# =============================================================================
# import seperate_star_gal
# seperate_star_gal.seperate(sextractorFile, img, starGal_out_file, plotLoc, psfMin, psfMax, fluxMin, fluxMax, psfMin_2, psfMax_2, fluxMin_2, fluxMax_2)
# =============================================================================

import coadd_detect
coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)

# =============================================================================
# import makeEandBmodes
# makeEandBmodes.EandB( '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy',  '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_r.npy', '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_i.npy', None, ebModeFile, img, '/scratch/bell/dutta26/wiyn_sim1/r_withMC.npy', '/scratch/bell/dutta26/wiyn_sim1/i_withMC.npy')
# 
# =============================================================================




# =============================================================================
# band ='i'
# img = '/scratch/bell/dutta26/wiyn_sim1/wted_coadds/'+band+'_coadd_wt.fits'
# med_img = '/scratch/bell/dutta26/wiyn_sim1/median_coadds/'+band+'_coadd_med_crop.fits'
# sextractorFile = '/scratch/bell/dutta26/wiyn_sim1/coadd_sim.cat'
# starGal_out_file = '/home/dutta26/codes/wiyn_wl_sim1/source_list.pk1'
# coadd_detect_file = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_' + str(band)+'.npy'
# ebModeFile = '/scratch/bell/dutta26/wiyn_sim1/master_arr_coadd.npy'
# plotLoc = '/scratch/bell/dutta26/wiyn_sim1/plot/phosim/'
# bandLoc = '/scratch/bell/dutta26/wiyn_sim1/'
# ir_band_data = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy'
# 
# import coadd_detect
# coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)
# 
# 
# band ='r'
# img = '/scratch/bell/dutta26/wiyn_sim1/wted_coadds/'+band+'_coadd_wt.fits'
# med_img = '/scratch/bell/dutta26/wiyn_sim1/median_coadds/'+band+'_coadd_med_crop.fits'
# sextractorFile = '/scratch/bell/dutta26/wiyn_sim1/coadd_sim.cat'
# starGal_out_file = '/home/dutta26/codes/wiyn_wl_sim1/source_list.pk1'
# coadd_detect_file = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_' + str(band)+'.npy'
# ebModeFile = '/scratch/bell/dutta26/wiyn_sim1/master_arr_coadd.npy'
# plotLoc = '/scratch/bell/dutta26/wiyn_sim1/plot/phosim/'
# bandLoc = '/scratch/bell/dutta26/wiyn_sim1/'
# ir_band_data = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy'
# 
# import coadd_detect
# coadd_detect.detect(starGal_out_file, img, ir_band_data, band, coadd_detect_file, plotLoc, bandLoc)
# 
# =============================================================================




