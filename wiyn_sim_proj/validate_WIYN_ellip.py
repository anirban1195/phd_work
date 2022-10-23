#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:49:11 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import helper 
import measure_pythonV
import matplotlib.pyplot as plt
from scipy.special import erf
import os,shutil

store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_chip10k50k.npy')
#store_long = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_long.npy')
coadd_file = '/scratch/halstead/d/dutta26/abell_2390/sample1.fits'
#os.remove('/scratch/halstead/d/dutta26/abell_2390/test1.fits')
#os.remove('/scratch/halstead/d/dutta26/abell_2390/test.fits')
#shutil.copy(coadd_file,'/scratch/halstead/d/dutta26/abell_2390/test1.fits')
#shutil.copy(coadd_file,'/scratch/halstead/d/dutta26/abell_2390/test.fits')

#f_test=fits.open('/scratch/halstead/d/dutta26/abell_2390/test.fits', mode='update')
min_dist_arr = np.load('/scratch/halstead/d/dutta26/abell_2390/min_dist.npy')
a,b = np.shape(store)
loc = np.where(((store[:,7]-store[:,38])<0) | ((store[:,8]-store[:,39])<0) )[0]
raArr=[]
decArr=[]
xArr=[]
yArr=[]
loc1=[]
count1Arr = []
count2Arr = [] 
count3Arr=[]
minLt = np.array([1, 20, 50 , 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000])
maxLt = np.array([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000,50000, 100000, 200000,500000, 1e6])

#minLt = np.array([2000])
#maxLt = np.array([9000])
#for k in range(len(minLt)):
for k in range(len(minLt)):
    count1 =0
    count2 = 0
    for idx in np.arange(a):
        if(store[idx,2] != 1  or store[idx,3]<= minLt[k] or  store[idx,3] > maxLt[k] or min_dist_arr[idx] < 20 ):
            continue
# =============================================================================
#         count1Arr.append(store[idx,7])
#         count2Arr.append(store[idx,7]-store[idx,38])
#         count3Arr.append(store[idx,3])
#     if(len(count1Arr)> 100):
#         n, bins, patches = plt.hist(x=count1Arr, bins=50, color='b',
#                                alpha=0.7, rwidth=0.85, density=True)
#         n, bins, patches = plt.hist(x=count2Arr, bins=50, color='r',
#                                alpha=0.7, rwidth=0.85, density=True)
#         plt.title("Flux range from "+str(minLt[k]) + ' - '+str(maxLt[k]) )
# =============================================================================
        
        if( ((store[idx,7]-store[idx,38])<0) | ((store[idx,8]-store[idx,39])<0)):
            count1 += 1
        else:
            count2 += 1
            
    count1Arr.append(count1)
    count2Arr.append(count2)
    print (count1, count2)
    
count1Arr =np.array(count1Arr)
count2Arr =np.array(count2Arr)

avgFlux = (maxLt+minLt)/2
plt.plot(np.log10(avgFlux), count1Arr/(count1Arr+count2Arr), 'g--')
# =============================================================================
#     if(np.isnan(store[idx,3]) or store[idx,3] == None or store[idx,2] == 1  or store[idx,3]<= 0.01 or min_dist_arr[idx] < 20 ):
#         continue
#     loc1.append(idx)
#     raArr.append(store[idx,0])
#     decArr.append(store[idx,1])
#     
#     x=int(store[idx,10])
#     y=int(store[idx,11])
#     f_test[0].data[y-20, x-20:x+20] = 1000
#     f_test[0].data[y+20, x-20:x+20] = 1000
#     f_test[0].data[y-20:y+20, x-20] = 1000
#     f_test[0].data[y-20:y+20, x+20] = 1000
# 
# f_test.flush()
# 
# =============================================================================


# =============================================================================
# e1_diff=[]
# e2_diff=[]
# 
# 
# for idx in loc1:
# # =============================================================================
# #     sigxx = store[idx,7]-store[idx,38]
# #     sigyy = store[idx,8]-store[idx,39]
# #     sigxy = store[idx,9]-store[idx,40]
# #     e1 = (sigxx-sigyy)/(sigxx+sigyy)
# #     e2 = 2*sigxy/(sigxx+sigyy)
# #     
# #     sigxx_l = store_long[idx,7]-store_long[idx,38]
# #     sigyy_l = store_long[idx,8]-store_long[idx,39]
# #     sigxy_l = store_long[idx,9]-store_long[idx,40]
# #     e1_l = (sigxx_l-sigyy_l)/(sigxx_l+sigyy_l)
# #     e2_l = 2*sigxy_l/(sigxx_l+sigyy_l)
# #     if(np.abs(e1) > 1 or np.abs(e2)>1 or np.abs(e1_l) > 1 or np.abs(e2_l)>1):
# #         continue
# #     if((e1-e1_l)>0.8 and store[idx,3]> 2000):
# #         print (store[idx,10], store[idx,11],store[idx,41],store[idx,42], e1,e1_l)
# # =============================================================================
#      
#     
#     size = np.sqrt(store[idx,7]+ store[idx,8])/np.pi
#     err_interp =0.35 
#     err_meas = np.sqrt((size / store[idx,3]) + ((4*size*size*np.pi*2000)/store[idx,3]**2))
#     err = np.sqrt(err_interp**2 + err_meas**2)     
#     print (err)             
#     sigxx = store[idx,7]-store[idx,38]+1
#     sigyy = store[idx,8]-store[idx,39]+1
#     sigxy = store[idx,9]-store[idx,40]
#     if(sigxx < 0):
#         A = store[idx,7]
#         B = store[idx,38]
#         s = err
#         sigxx = 0.5*(A-B)*(erf((A-B)/np.sqrt(2)/s) + 1)   + s/np.sqrt(2*np.pi)*np.exp(-(A-B)**2/2/s/s)
#     if(sigyy < 0):
#         A = store[idx,8]
#         B = store[idx,39]
#         s = err
#         sigyy = 0.5*(A-B)*(erf((A-B)/np.sqrt(2)/s) + 1)   + s/np.sqrt(2*np.pi)*np.exp(-(A-B)**2 /2/s/s)
#     #C = 0.5*(A-B)*(erf((A-B)/sqrt(2)/s) + 1)   + s/sqrt(2*PI)*exp(-(A-B)^2/2/s/s)
# 
#             
#     e1 = (sigxx-sigyy)/(sigxx+sigyy)
#     e2 = 2*sigxy/(sigxx+sigyy)
#     
#     sigxx_l = store_long[idx,7]-store_long[idx,38]
#     sigyy_l = store_long[idx,8]-store_long[idx,39]
#     sigxy_l = store_long[idx,9]-store_long[idx,40]
#     e1_l = (sigxx_l-sigyy_l)/(sigxx_l+sigyy_l)
#     e2_l = 2*sigxy_l/(sigxx_l+sigyy_l)
#     if(np.abs(e1) > 1 or np.abs(e2)>1 or np.abs(e1_l) > 1 or np.abs(e2_l)>1):
#         continue
#     #if((e1-e1_l)>0.8 and store[idx,3]> 2000):
#     #    print (store[idx,10], store[idx,11],store[idx,41],store[idx,42], e1,e1_l)    
#         
#         
#     
#     e1_diff.append(e1-e1_l)
#     e2_diff.append(e2-e2_l)
#     #plt.plot(store[idx,3],np.abs(e1-e1_l), 'b.' )
# n, bins, patches = plt.hist(x=e1_diff, bins=20, color='b',
#                            alpha=0.5, rwidth=0.85, density=True)   
# n, bins, patches = plt.hist(x=e2_diff, bins=20, color='r',
#                            alpha=0.5, rwidth=0.85, density=True) 
# 
# 
# 
# =============================================================================
