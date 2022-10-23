#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 08:31:00 2021

@author: dutta26
"""


from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
import helper
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess

filename = '/scratch/halstead/d/dutta26/abell_2390_limited/abell_ir_coadd_wted1.fits'
#filename = str(sys.argv[1])
#Read the catalog 
catalog = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test.cat'
#catalog = str(sys.argv[2])
with open(catalog) as f:
    content = f.readlines()
    
Raindex=5
Decindex=6
Aindex = 7
Bindex = 8
thetaindex = 9


#Make an array containing x and y indices 
raList =[]
decList =[]
sizeList=[]
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    raList.append(float(content[j].split()[Raindex])) 
    decList.append( float(content[j].split()[Decindex])) 
    A=float(content[j].split()[Aindex])
    theta = float(content[j].split()[thetaindex])  
    sizeList.append(np.abs(3*A*np.cos(theta)))
    

#Falg is set to 1 if any nearby objects
nearFlagList =[]

xList, yList = helper.convertToXY(raList, decList, filename)
# =============================================================================
# #Conversion to x and y complete
# #Setting nearby flags
# for j in range(len(xList)):
#     
#     lowLimit = j-1000
#     if(lowLimit<0):
#         lowLimit = 0
#     highLimit = j + 1000
#     if(highLimit> len(xList)):
#         highLimit = len(xList)
#     nearbyFlag = 0
#     
#     for k in range(lowLimit, highLimit):
#         if(j == k):
#             continue
#         x1 = xList[j]
#         x2 = xList[k]
#         y1 = yList[j]
#         y2 = yList[k]
#         totSize = sizeList[j]+ sizeList[k]
#         sepn = np.sqrt((x1-x2)**2 + (y1-y2)*2)
#         if(totSize >= sepn):
#             nearbyFlag = 1
#     nearFlagList.append(nearbyFlag)
#             
# =============================================================================
# # =============================================================================
# fileList =[]
# loc = '/scratch/halstead/d/dutta26/abell_2390_limited/r/'
# for fitsFiles in os.listdir(loc):
#     fileList.append(loc+fitsFiles)
# #fileList = '/scratch/halstead/d/dutta26/abell_2390_limited/r/'+os.listdir('/scratch/halstead/d/dutta26/abell_2390_limited/r/')  
# swarpLoc = '/home/dutta26/apps/bin/bin/'
# 
# for files in fileList:
#     #Run sextractor 
#     tempOut = '/scratch/halstead/d/dutta26/abell_2390_limited/coadd_temp.fits'
#     tempOut = '/scratch/halstead/d/dutta26/abell_2390_limited/abell_ir_coadd_wted1.fits'
#     tempWtOut = '/scratch/halstead/d/dutta26/abell_2390_limited/coadd_temp.weight.fits'
# # =============================================================================
# #     bashCommand = './swarp '+files + ' -c /home/dutta26/default_psfMap.swarp -IMAGEOUT_NAME ' +tempOut + ' -WEIGHTOUT_NAME '+tempWtOut
# #     process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
# #     output, error = process.communicate()
# # =============================================================================
#     
#     xList, yList = helper.convertToXY(raList, decList, tempOut)
#     indexList = helper.detectBad(tempOut, xList, yList, sizeList)
#     
#     fluxList, psfList, bkgList, e1List, e2List, muXList, muYList = helper.runMeasure(tempOut, xList, yList, indexList, sizeList)
#     
#     break
#     
#     
# =============================================================================
c=0    
for j in range(len(xList)):
    for k in range(len(yList)):
        if (xList[j] >5940 and xList[j] <5960 and yList[k] > 11688 and yList[k]<11708):
            c=c+1 
            
    

