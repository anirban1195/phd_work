#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:29:05 2023

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt

master_frame = np.load('/scratch/bell/dutta26/wiyn_sim/master_arr_coaddMC.npy')
ir_coadd_data = np.load('/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy')


size = np.sqrt(ir_coadd_data[:,7] + ir_coadd_data[:,8])
loc = np.where(  (master_frame[:,4]> 1) & (ir_coadd_data[:,2]==0) & (ir_coadd_data[:,3]>10)
               & (ir_coadd_data[:,82]==0))[0]

errArr1=[]
errArr2=[]
errArr3=[]
eArr1=[]
eArr2=[]
sizeArr1=[]
sizeArr2=[]
cnt =0
for ind in loc:
    size = np.sqrt(master_frame[ind,6] + master_frame[ind,7])
    e1 = master_frame[ind,2]
    e2 = master_frame[ind,3]
# =============================================================================
#     if(abs(e1)>1 or abs(e2)>1   ):
#         continue
# =============================================================================
    
    if(size>0.15 and size<2.5):
        print(ind)
        errArr1.append(master_frame[ind, 10])
        eArr1.append(e2**2 + e1**2)
        sizeArr1.append(size)
    if(size>2.5 and size<80):
        errArr2.append(master_frame[ind, 10])
        eArr2.append(e2**2 + e1**2)
        sizeArr2.append(size)
# =============================================================================
#     if(size>5 and size<8):
#         errArr3.append(master_frame[ind, 15])
# =============================================================================
        
# =============================================================================
# n, bins, patches = plt.hist(x=errArr1, bins='auto',histtype=u'step', color='r', label='Small', density = True)
# n, bins, patches = plt.hist(x=errArr2, bins='auto',histtype=u'step', color='b', label='Medium', density = True)
# n, bins, patches = plt.hist(x=errArr3, bins='auto',histtype=u'step', color='k', label='Large', density = True)
# plt.legend()
# =============================================================================

print (len(errArr1)+ len(errArr2)+len(errArr3))
#plt.plot(eArr1, errArr1, 'r+', markersize= 2, label = 'Small Source')
plt.plot(eArr2, errArr2, 'b+', markersize= 2, label = 'Medium Source')
plt.xlabel('Ellipticity')
plt.ylabel('Error in Elliptcity^2')
plt.plot(np.arange(0,100,0.1), np.arange(0,100,0.1), 'k--')
plt.legend()