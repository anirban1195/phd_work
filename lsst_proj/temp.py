#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 10:46:37 2022

@author: dutta26
"""
import os,shutil,sys
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import helper, helper1, measure_pythonV,sys
# =============================================================================
# for folders in os.listdir('/scratch/halstead/d/dutta26/lsst/3/'):
#     for files in os.listdir('/scratch/halstead/d/dutta26/lsst/3/'+folders):
#         if('CELL' in files):
#             continue
#         else:
#             shutil.copy('/scratch/halstead/d/dutta26/lsst/3/'+folders+'/'+files, '/scratch/halstead/d/dutta26/lsst/work3/'+files)
#             
#             
# =============================================================================
            
f=open('/home/dutta26/temp.ascii', 'w+')
# =============================================================================
# for files in os.listdir('/scratch/halstead/d/dutta26/lsst/work1/'):
#     if('weight' in files):
#         continue
#     f.write('/scratch/halstead/d/dutta26/lsst/work1/'+files+ '\n')
# =============================================================================
# =============================================================================
# for files in os.listdir('/scratch/halstead/d/dutta26/lsst/work2/'):
#     if('weight' in files):
#         continue
#     f.write('/scratch/halstead/d/dutta26/lsst/work2/'+files+ '\n')
# =============================================================================
for files in os.listdir('/scratch/halstead/d/dutta26/lsst/work3/'):
    if('weight' in files):
        continue
    f.write('/scratch/halstead/d/dutta26/lsst/work3/'+files+ '\n')
f.close()
# =============================================================================
# import numpy as np 
# import helper
# xList = np.arange(2470, 4075, 177.66)
# yList = np.arange(2207, 3759, 171.5)
# raArr=[]
# decArr=[]
# xArr=[]
# yArr=[]
# for j in range(len(xList)):
#     for k in range(len(yList)):
#         yArr.append(yList[j])
#         xArr.append(xList[k])
#         
# raArr,decArr = helper.convertToRaDec(xArr, yArr, '/scratch/halstead/d/dutta26/lsst/star_grid_bkg.fits')
#         
# np.save('/scratch/halstead/d/dutta26/lsst/raArr_wiyn.npy', raArr)
# np.save('/scratch/halstead/d/dutta26/lsst/decArr_wiyn.npy', decArr)
# =============================================================================

# =============================================================================
# allFiles = os.listdir('/scratch/halstead/d/dutta26/lsst/work1/')
# for files in allFiles:
#     if('weight' in files):
#         continue
#     index =''
#     for j in range(14, 16):
#         if(files[j] == '_'):
#             break
#         else:
#             index += files[j]
#     
#     print (index)
#     catFile = '/scratch/halstead/d/dutta26/lsst/wiyn_cat1/'+index+'.txt1_flat'
#     f=open(catFile)
#     content = f.readlines()
#     f.close()
#     see = float(content[6][7:12])
#     
#     #Find ZP
#     raList = np.load ('/scratch/halstead/d/dutta26/lsst/raArr_wiyn.npy')
#     decList = np.load ('/scratch/halstead/d/dutta26/lsst/decArr_wiyn.npy')  
#     
#     
#     #Copy to a weight files
#     wtFile = '/scratch/halstead/d/dutta26/lsst/work1/'+files[:-5]+'.weight.fits'
#     shutil.copy('/scratch/halstead/d/dutta26/lsst/work1/'+files, wtFile)
#     f=fits.open(wtFile, mode = 'update')
#     data = f[0].data
#     xList, yList = helper.convertToXY(raList, decList, wtFile)
#     for j in range(96,97):
#         y=int(yList[j])
#         x=int(xList[j])
#         cut = data[y-20:y+20, x-20:x+20]
#         flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#         break
#     zp = 16.2 + 2.5*np.log10(flux/100)
#     
#     bkg_std = sigma_clipped_stats(data)[2]
#     f[0].data = np.ones(np.shape(data)) * (100* 10**(zp-27)) /(bkg_std**2 * see**2) 
#     f.flush()
#     print (files, see, zp,(100* 10**(zp-27)) /(bkg_std**2 * see**2))
#     #sys.exit()
# 
# allFiles = os.listdir('/scratch/halstead/d/dutta26/lsst/work2/')    
# for files in allFiles:
#     if('weight' in files):
#         continue
#     index =''
#     for j in range(14, 16):
#         if(files[j] == '_'):
#             break
#         else:
#             index += files[j]
#     
#     print (index)
#     catFile = '/scratch/halstead/d/dutta26/lsst/wiyn_cat1/'+index+'.txt2_flat'
#     f=open(catFile)
#     content = f.readlines()
#     f.close()
#     see = float(content[6][7:12])
#     
#     #Find ZP
#     raList = np.load ('/scratch/halstead/d/dutta26/lsst/raArr_wiyn.npy')
#     decList = np.load ('/scratch/halstead/d/dutta26/lsst/decArr_wiyn.npy')  
#     
#     
#     #Copy to a weight files
#     wtFile = '/scratch/halstead/d/dutta26/lsst/work2/'+files[:-5]+'.weight.fits'
#     shutil.copy('/scratch/halstead/d/dutta26/lsst/work2/'+files, wtFile)
#     f=fits.open(wtFile, mode = 'update')
#     data = f[0].data
#     xList, yList = helper.convertToXY(raList, decList, wtFile)
#     for j in range(96,97):
#         y=int(yList[j])
#         x=int(xList[j])
#         cut = data[y-20:y+20, x-20:x+20]
#         flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#         break
#     zp = 16.2 + 2.5*np.log10(flux/100)
#     
#     bkg_std = sigma_clipped_stats(data)[2]
#     f[0].data = np.ones(np.shape(data)) * (100* 10**(zp-27)) /(bkg_std**2 * see**2) 
#     f.flush()
#     print (files,see,zp, (100* 10**(zp-27)) /(bkg_std**2 * see**2))
#     
# allFiles = os.listdir('/scratch/halstead/d/dutta26/lsst/work3/')    
# for files in allFiles:
#     if('weight' in files):
#         continue
#     index =''
#     for j in range(14, 16):
#         if(files[j] == '_'):
#             break
#         else:
#             index += files[j]
#     
#     print (index)
#     catFile = '/scratch/halstead/d/dutta26/lsst/wiyn_cat1/'+index+'.txt3_flat'
#     f=open(catFile)
#     content = f.readlines()
#     f.close()
#     see = float(content[6][7:12])
#     
#     #Find ZP
#     raList = np.load ('/scratch/halstead/d/dutta26/lsst/raArr_wiyn.npy')
#     decList = np.load ('/scratch/halstead/d/dutta26/lsst/decArr_wiyn.npy')  
#     
#     
#     #Copy to a weight files
#     wtFile = '/scratch/halstead/d/dutta26/lsst/work3/'+files[:-5]+'.weight.fits'
#     shutil.copy('/scratch/halstead/d/dutta26/lsst/work3/'+files, wtFile)
#     f=fits.open(wtFile, mode = 'update')
#     data = f[0].data
#     xList, yList = helper.convertToXY(raList, decList, wtFile)
#     for j in range(97,98):
#         y=int(yList[j])
#         x=int(xList[j])
#         cut = data[y-28:y+28, x-28:x+28]
#         flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#         break
#     zp = 16.2 + 2.5*np.log10(flux/100)
#     
#     bkg_std = sigma_clipped_stats(data)[2]
#     f[0].data = np.ones(np.shape(data)) * (100* 10**(zp-27)) /(bkg_std**2 * see**2) 
#     f.flush()
#     print (files, see, zp, (100* 10**(zp-27)) /(bkg_std**2 * see**2))
# 
# 
# =============================================================================
