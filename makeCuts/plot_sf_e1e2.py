#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 08:10:41 2023

@author: dutta26
"""

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import sys
r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r.pk1'))
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')

#For r band do flux 50-60, 10-12, 1-1.25

loc = np.where((r_coadd_df[:,3]> 40) & (r_coadd_df[:,3]<50) & (r_coadd_df[:,2]!= 1))[0]
fig, axs = plt.subplots(2)
fig.suptitle('Flux = 1-1.2. SNR in frames ~1 ')
a,b,c = np.shape(r_sf_df)
r = np.zeros((10))
s = np.zeros((16))
testArr =[]

for i in loc[20:300]:
    xx_coadd = r_coadd_df[i, 35] - r_coadd_df[i, 38]
    yy_coadd =  r_coadd_df[i, 36] - r_coadd_df[i, 39]
    xy_coadd =  r_coadd_df[i, 37] - r_coadd_df[i, 40]
    e1_coadd = (xx_coadd - yy_coadd)/(xx_coadd + yy_coadd)
    e2_coadd = 2*xy_coadd / (xx_coadd+yy_coadd)
    
    e1Arr=[]
    e2Arr=[]
    for j in range(a-1):
        if(r_sf_df[j,i, 38]< 0 or r_sf_df[j,i,39]< 0):
            r[1] += 1
            continue
        if(r_sf_df[j,i, 13]== 1 or r_sf_df[j,i,14] == 1 ):
            r[2] += 1
            
            continue
        if(1 in r_sf_df[j,i, 54:71]):
            s += r_sf_df[j,i, 54:71]
            if(r_sf_df[j,i, 62] == 1):
                print (j,i)
            r[3] += 1
            continue
        if(r_sf_df[j,i, 47]<0 or r_sf_df[j,i, 48]<0):
            r[4] += 1
            continue
        if(r_sf_df[j,i, 51] == 0):
            r[5] += 1
            sys.exit()
        xx = r_sf_df[j,i, 51] - r_sf_df[j,i, 38]
        yy = r_sf_df[j,i, 52] - r_sf_df[j,i, 39]
        xy = r_sf_df[j,i, 53] - r_sf_df[j,i, 40]
        e1 = (xx - yy)/ (xx+yy)
        e2 = 2*xy/(xx+yy)
        #print (r_sf_df[j,i, 51], r_sf_df[j,i, 38], j)
        e1Arr.append(e1)
        e2Arr.append(e2)
        
    
    #print (len(e1Arr))
    testArr.append(len(e1Arr))
    
    
    axs[0].plot( e1_coadd,  np.nanmedian(e1Arr), 'b.', markersize = 3)
    axs[1].plot(e2_coadd, np.nanmedian(e2Arr),  'b.', markersize = 3)
#axs[0].xlabel('e1 Median of Single Frames')
axs[0].set(xlabel='e1_coadd ', ylabel='e1 Median of Single Frames')
axs[1].set(xlabel='e2_coadd ', ylabel='e2 Median of Single Frames')
print (np.median(testArr)   )     
print (r/((a-1)*len(loc)) )
print (s/((a-1)*len(loc)) )