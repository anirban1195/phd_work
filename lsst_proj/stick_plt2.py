#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 08:03:38 2022

@author: dutta26
"""


import matplotlib.pyplot as plt


import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import griddata
import sys,helper,os

band = '5'
#ir_coadd_data = pd.read_pickle('/scratch/halstead/d/dutta26/m_38/df_'+band+'.pk1')
#store = np.array(ir_coadd_data)

store = np.load('/scratch/halstead/d/dutta26/m_38/singleFrame_'+band+'.npy')
imp_xx =[]
imp_yy=[]
imp_xy =[]

med_xx=[]
med_yy=[]
med_xy=[]

for k in range(2,20):
    loc = np.where((store[k,:,2] == 1) & (store[k,:,3] > 40000) &  (store[k,:,3] < 50000000) 
                   & (store[k,:,12] == 99) &(store[k,:,14] == 0))[0]
    # =============================================================================
    # ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
    # store = np.array(ir_coadd_data)
    # loc = np.where((store[:,2] == 1) & (store[:,3] > 500) & (store[:,7] < 8) &(store[:,7] > 3))[0]
    # # =============================================================================
    # =============================================================================
    # store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
    # loc= (np.where((store[:,2] == 1) & (store[:,3] > 2000)  &(store[:,3] < 40000) 
    #                & (store[:,13] == 1) & (store[:,7] >3) & (store[:,7] <6) ))[0]
    # =============================================================================
    a=store[k,loc,7]
    b=store[k,loc,7]- store[k,loc,38]
    xList=store[k,loc,10]
    yList=store[k,loc, 11]
    sizeList = np.sqrt(store[k,loc,7] + store[k,loc,8])
    
    sigxx = sigma_clipped_stats(store[k,loc,7])[2]
    sigxx_imp = sigma_clipped_stats(store[k,loc,7]-  store[k,loc,38])[2]
    #print (sigxx, sigxx_imp,(sigxx - sigxx_imp)/sigxx)
    imp_xx.append((sigxx - sigxx_imp)/sigxx)
    med_xx.append(sigma_clipped_stats(store[k,loc,7])[2])
    
    sigyy = sigma_clipped_stats(store[k,loc,8])[2]
    sigyy_imp = sigma_clipped_stats(store[k,loc,8]-  store[k,loc,39])[2]
    #print (sigyy, sigyy_imp,(sigyy - sigyy_imp)/sigyy)
    imp_yy.append((sigyy - sigyy_imp)/sigyy)
    med_yy.append(sigma_clipped_stats(store[k,loc,8])[2])
    
    sigxy = sigma_clipped_stats(store[k,loc,9])[2]
    sigxy_imp = sigma_clipped_stats(store[k,loc,9]-  store[k,loc,40])[2]
    imp_xy.append((sigxy - sigxy_imp)/sigxy)
    med_xy.append(sigma_clipped_stats(store[k,loc,9])[2])
    
    print (sigma_clipped_stats(store[k,loc,8]))
    print (sigma_clipped_stats(store[k,loc,8]-  store[k,loc,39]))
    print ('***********')
# =============================================================================
#     #e1 = (store[k,loc,7] - store[k,loc,8])/(store[k,loc,7] + store[k,loc,8])
#     #e2 = 2* store[k,loc,9] /(store[k,loc,7] + store[k,loc,8])
#     sigmaxx = store[k,loc,7] 
#     sigmayy = store[k,loc,8] 
#     sigmaxy = store[k,loc,9] 
#     
#     print (sigma_clipped_stats(sigmaxx))
#     print (sigma_clipped_stats(sigmayy))
#     print (sigma_clipped_stats(sigmaxy))
#     
#     e1 = (sigmaxx-sigmayy)/(sigmaxx+sigmayy)
#     e2 = 2*sigmaxy / (sigmaxx+sigmayy)
#     
#     print (sigma_clipped_stats(e1))
#     print (sigma_clipped_stats(e2))
#     
#     ellip = np.sqrt(e1**2+e2**2)
#     theta = 0.5*np.arctan2(e2,e1)
#     cosArr = np.cos(theta)
#     sinArr = np.sin(theta)
#     
#     fact=5
#     #Get the rage of values of x and y
#     xShape = int(17000/fact) + 10
#     yShape = int(17000/fact) + 10
#     xList= xList/fact
#     yList = yList/fact
#     
#     #f=fits.open('/scratch/halstead/d/dutta26/m_38/temp.fits', mode='update')
#     #data = np.array(f[0].data)
#     #img = np.zeros((yShape, xShape), dtype=np.float32)
#     count = 0
#     for j in range(len(xList) - 1):
#         dist = (ellip[j]/ 0.01) *15
#         if(dist > 100):
#             continue
#         xNew = xList[j] + dist*cosArr[j]
#         yNew = yList[j] + dist*sinArr[j]
#         xArr = [int(xList[j]), int(xNew)]
#         yArr = [int(yList[j]), int(yNew)]
#         plt.plot(xArr, yArr, 'b')
#         count += 1
#         #data[int(fact*yList[j] - 20):int(fact*yList[j] + 20), int(fact*xList[j] - 20)] = 100000
#         #data[int(fact*yList[j] - 20):int(fact*yList[j] + 20), int(fact*xList[j] + 20)] = 100000
#         #data[int(fact*yList[j] - 20), int(fact*xList[j] - 20):int(fact*xList[j] + 20)] = 100000
#         #data[int(fact*yList[j] + 20), int(fact*xList[j] - 20):int(fact*xList[j] + 20)] = 100000
#     #f[0].data = data
#     #f.flush()
#     print (count)
#     sys.exit()
#     
# =============================================================================
    
# =============================================================================
#     points=[]
#     values =[] 
#     grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
#     
#     for j in range(len(xList)- 1):
#         x = int(round(xList[j]))
#         y = int(round(yList[j]))
#         
#         points.append([y,x])
#         values.append( e1[j])
#     
#     
#     grid_z1 = griddata(points, values, (grid_y, grid_x), method='linear')    
#     grid_z1 = np.array(grid_z1, dtype = np.float32).T
#     hdu = fits.PrimaryHDU(grid_z1) 
#     hdu.writeto('/scratch/halstead/d/dutta26/m_38/filter'+band+'/data/ellip_'+str(k)+'.fits', clobber=True)     
# 
# =============================================================================
