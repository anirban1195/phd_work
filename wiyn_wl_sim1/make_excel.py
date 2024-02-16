#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 08:26:01 2023

@author: dutta26
"""
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os,sys
import pandas as pd
import helper, helper1
from astropy.stats import sigma_clipped_stats
import subprocess
band_coadd_npy = np.load('/home/dutta26/codes/wiyn_wl_sim/coaddSc_r.npy')
band_sf_npy = np.load('/scratch/bell/dutta26/wiyn_sim/r_withMC.npy')

loc = np.where((band_sf_npy[0,:,12] == 99) & (band_sf_npy[0,:,3] >0) & (band_sf_npy[0,:,2] ==1))
indList= [1764, 3363, 17396]
#indList= [1764, 3363, 1387]

columns_coadd = ['ra', 'dec', 'star_flag',
'flux', 'mux', 'muy', 'bkg', 
'xx', 'yy', 'xy','x','y','force_flag', 
 'vert_flag',  'bad_flag',  'back_sky',  'seeing',  'zp',  'fwhm',  'mjd' , 'airmass', 
 'mphase',  'mAngle' , 'expTime' , 'focus' , 'zp_n' , 'skymag' , 'depth' , 'mRa' ,
 'mDec', 'scale_fact', 'flux_sf', 'mux_sf', 'muy_sf', 'bkg_sf', 'xx_sf', 'yy_sf', 'xy_sf',
 'interp_xx', 'interp_yy', 'interp_xy', 'std_xx', 'std_yy', 'std_xy','file_id',
 'bad_loc', 'bad_loc_1', 'ccd_n0', 'chip_x',  'chip_y''Empty','sig_xx','sig_yy','sig_xy','flag1 ',
'flag2','flag3','flag4','Empty','Empty','flag5','flag6','flag7','flag8','flag9',
'flag10','flag11','flag12','e_xx','e_yy','e_xy ','e_xx','e_yy','e_xy ','Error',
'Error','Error','success no.','boundary radius','overlap flag (outdated)','closest centroid ','angle closest',
'bkg flag','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty',
'Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty']



store0= np.zeros((31,100), dtype= np.float32)
store1= np.zeros((31,100), dtype= np.float32)
store2= np.zeros((31,100), dtype= np.float32)

store0[0,:] = band_coadd_npy[indList[0],:]
store1[0,:] = band_coadd_npy[indList[1],:]
store2[0,:] = band_coadd_npy[indList[2],:]

store1[0,35] -= store1[0,38]
store1[0,36] -= store1[0,39]
store1[0,37] -= store1[0,40]

for j in range(30):
    store0[j+1,0:77] = band_sf_npy[j, indList[0], :]
    store0[j+1,78] = band_sf_npy[j, indList[0], 31]*band_sf_npy[j, indList[0], 30]
    store0[j+1,79] = band_sf_npy[j, indList[0], 31]*band_sf_npy[j, indList[0], 30]
    
    store1[j+1,0:77] = band_sf_npy[j, indList[1], :]
    store1[j+1,78] = band_sf_npy[j, indList[1], 31]*band_sf_npy[j, indList[1], 30]
    
    
    store2[j+1,0:77] = band_sf_npy[j, indList[2], :]
    store2[j+1,78] = band_sf_npy[j, indList[2], 31]*band_sf_npy[j, indList[2], 30]
    
   

# =============================================================================
# df = pd.DataFrame(store, columns = columns_coadd)
# df.to_excel("/scratch/bell/dutta26/wiyn_sim/converged stars.xlsx", sheet_name='Sheet_name_'+str(count)  )  
# count += 1
# =============================================================================
df0 = pd.DataFrame(store0, columns = columns_coadd)
df1 = pd.DataFrame(store1, columns = columns_coadd)
df2 = pd.DataFrame(store2, columns = columns_coadd)
with pd.ExcelWriter("/scratch/bell/dutta26/wiyn_sim/forced_gal2.xlsx") as writer:  
    df0.to_excel(writer, sheet_name='Sheet_name_0')
    df1.to_excel(writer, sheet_name='Sheet_name_1')
    df2.to_excel(writer, sheet_name='Sheet_name_2')

