#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 16:32:39 2022

@author: dutta26
"""

from astropy.io import fits
import numpy as np 
import os,sys
import subprocess
import gzip
import matplotlib.pyplot as plt
import shutil
from astropy.stats import sigma_clipped_stats
import helper
import measure_pythonV
import pandas as pd
import matplotlib.pyplot as plt

f=open('/scratch/halstead/d/dutta26/lsst/lsst_catalog1/lsst.cat')
cont = f.readlines()
f.close()

raList = []
decList = []
magArr =[]
for j in range(len(cont)):
    if('star' not in cont[j]):
        continue
    a = cont[j].split()
    if(float(a[4])> 25):
        continue
    raList.append(float(a[2]))
    decList.append(float(a[3]))
    magArr.append(float(a[4]))
    
coadd_file = '/scratch/halstead/d/dutta26/lsst/coadd1_ir.fits'   
singleFrame = '/scratch/halstead/d/dutta26/lsst/filter4/1001/final_coadd.fits' 

xList, yList = helper.convertToXY(raList, decList, coadd_file)
xList1, yList1 = helper.convertToXY(raList, decList, singleFrame)
f=fits.open(coadd_file)
coadd_data = f[0].data
f.close()

f=fits.open(singleFrame)
sf_data = f[0].data
f.close()




ySize,xSize = np.shape(coadd_data)
ySize1,xSize1 = np.shape(sf_data)

fluxArr =[]
sizeArr =[]
sizeArr1 =[]
for j in range(len(xList1)):
    
    if(xList[j]>xSize-25 or xList[j]<25 or yList[j]>ySize-25 or yList[j]<25):
        continue
    if(xList1[j]>xSize1-25 or xList1[j]<25 or yList1[j]>ySize1-25 or yList1[j]<25):
        continue
    
    x=int(xList[j])
    y=int(yList[j])
    cut = coadd_data[y-15: y+15, x-15: x+15]
    
    x=int(xList1[j])
    y=int(yList1[j])
    cut1 = sf_data[y-15: y+15, x-15: x+15]

    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
    
    flux1, mux1, muy1, e1_1, e2_1, bkg1, psf1, sigxx1, sigyy1, sigxy1 = measure_pythonV.measure_v2(cut1)
    
    if(flux == None or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
        continue
    
    if(flux1 == None or e1_1 == None or e2_1 ==None or  np.isnan(e1_1) or np.isnan(e2_1) or np.isnan(psf1)):
        continue
    calc_flux = 10**((31 - magArr[j])/2.5)
    dev = (flux-calc_flux)/calc_flux
    dev1 = (flux1-calc_flux)/calc_flux
    if(abs(dev1)> 0.2 or abs(dev)>0.2):
        continue
    fluxArr.append(flux)
    sizeArr.append(psf)
    sizeArr1.append(psf1)

fluxArr = np.array(fluxArr)
sizeArr = np.array(sizeArr)
sizeArr1 = np.array(sizeArr1)

fluxUp = [100,200,400, 800, 2000, 5000, 10000, 40000, 100000, 200000, 400000]
fluxDown = [50,100, 200, 400, 800, 2000, 5000, 10000, 40000, 100000, 200000]
meanSize = sigma_clipped_stats(sizeArr)[0]
for j in range(len(fluxUp)):
    loc = np.where( (fluxArr< fluxUp[j]) & (fluxArr> fluxDown[j]))[0]
    print (len(loc))
    a =sigma_clipped_stats(sizeArr[loc])
    plt.errorbar(np.log10(fluxDown[j]), a[1], yerr= 0, fmt='.g')
    a =sigma_clipped_stats(sizeArr1[loc])
    plt.errorbar(np.log10(fluxDown[j]), a[1],yerr= 0, fmt='+g')
#plt.errorbar(0, 0, yerr= 0, fmt='+g', label = 'LSST Simulation Co-add')


