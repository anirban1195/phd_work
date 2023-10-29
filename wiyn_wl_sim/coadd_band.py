#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 10:57:55 2023

@author: dutta26
"""

import numpy as np
from astropy.io import fits
import helper_phosim,correct_wcs
import subprocess,sys

phosimLoc = '/home/dutta26/Downloads/phosim_release/'
filtArr =['u', 'g', 'r', 'i', 'z']
for filt in filtArr:
# =============================================================================
#     filt = 'ir'
#     idNo = 1000
#     noImages = 30
#     f=open('/scratch/bell/dutta26/wiyn_sim/coadd_cat.ascii', 'w+')
#     for j in range(noImages):
#         if(filt == 'u' and idNo == 1008): #Defective image
#             idNo += 1
#             continue
#         print ('*********************************************')
#         print ('/scratch/bell/dutta26/wiyn_sim/'+filt+'/'+str(idNo)+'/output.fits \n')
#         print ('********************************************started*')
#         #correct_wcs.correct_wcs(str(idNo), filt)
#         #helper_phosim.makeWeight(str(idNo), filt)
#         
#         
#         f.write('/scratch/bell/dutta26/wiyn_sim/'+filt+'/'+str(idNo)+'/output.fits \n')
#         idNo += 1
#         
#     f.close()
#     
#     #For weighted coadd
#     swarpLoc = '/home/dutta26/apps/bin/bin/'
#     swarpCommand = './swarp @/scratch/bell/dutta26/wiyn_sim/coadd_cat.ascii -c /scratch/bell/dutta26/wiyn_sim/default.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/wiyn_sim/wted_coadds/'+filt+'_coadd_wt.fits -WEIGHTOUT_NAME  /scratch/bell/dutta26/wiyn_sim/wted_coadds/'+filt+'_coadd_wt.weight.fits -COMBINE_TYPE WEIGHTED'
#     process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
#     output, error = process.communicate()
#         
# =============================================================================
    
# =============================================================================
#     #For Median coadd
#     swarpLoc = '/home/dutta26/apps/bin/bin/'
#     swarpCommand = './swarp @/scratch/bell/dutta26/wiyn_sim/coadd_cat.ascii -c /scratch/bell/dutta26/wiyn_sim/default.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/wiyn_sim/median_coadds/'+filt+'_coadd_med.fits -WEIGHTOUT_NAME  /scratch/bell/dutta26/wiyn_sim/median_coadds/'+filt+'_coadd_med.weight.fits -COMBINE_TYPE MEDIAN'
#     
#     process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
#     output, error = process.communicate()
# =============================================================================

    #For the NFW profiles
    loc = '/scratch/bell/dutta26/backup/'
    idNo = 1000
    noImages = 30
    f=open('/scratch/bell/dutta26/wiyn_sim/coadd_cat.ascii', 'w+')
    for j in range(noImages):
        if(filt == 'u' and idNo == 1008): #Defective image
            idNo += 1
            continue
        print ('*********************************************')
        print ('scratch/bell/dutta26/backup/'+filt+'/'+str(idNo)+'/output.fits \n')
        print ('********************************************started*')
        #correct_wcs.correct_wcs(str(idNo), filt, loc)
        #helper_phosim.makeWeight(str(idNo), filt, loc)
        
        
        f.write('/scratch/bell/dutta26/backup/'+filt+'/'+str(idNo)+'/output.fits \n')
        idNo += 1
        
    f.close()
    
    #For weighted coadd
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    swarpCommand = './swarp @/scratch/bell/dutta26/wiyn_sim/coadd_cat.ascii -c /scratch/bell/dutta26/wiyn_sim/default.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/backup/wted_coadds/'+filt+'_coadd_wt.fits -WEIGHTOUT_NAME  /scratch/bell/dutta26/backup/wted_coadds/'+filt+'_coadd_wt.weight.fits -COMBINE_TYPE WEIGHTED'
    process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
        
        
    