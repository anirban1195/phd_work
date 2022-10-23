#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 09:09:37 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import wquantiles,sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime

ir_coadd_data = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_ir.pk1')
ir_coadd_data = np.array(ir_coadd_data)
zArr = np.load('/scratch/halstead/d/dutta26/lsst/lsst_zArr.npy')


#galLoc =  np.where(ir_coadd_data[:,41] > 0)[0]
master_frame = np.zeros((len(ir_coadd_data), 5), dtype =np.float32)
#Converge the measurements of 5 frames to 1 single frame 


master_frame[:,0] = ir_coadd_data[:,7]
master_frame[:,1] = ir_coadd_data[:,8]
master_frame[:,2] = ir_coadd_data[:,9]
master_frame[:,3] = ir_coadd_data[:, 10]
master_frame[:,4] = ir_coadd_data[:, 11]
    

       
coadd_file = '/scratch/halstead/d/dutta26/lsst/coadd1_ir.fits'
f=fits.open(coadd_file)
data = np.array(f[0].data)
ySize,xSize = np.shape(data)
f.close()  
del data 
    
chopSize = 50
alphax = 500

mod_ySize = int(round(ySize/chopSize)) + 2
mod_xSize = int(round(xSize/chopSize)) + 2
imgE = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
imgB = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)


#master_i = master_frame[0,:,:]
#master_r = master_frame[1,:,:]
cd_1 = np.array(pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_1.pk1'))
sf_1 = np.load('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_1.npy')

cd_2 =np.array(pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_2.pk1'))
sf_2 =np.load('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_2.npy')

cd_3 =np.array(pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_3.pk1'))
sf_3 =np.load('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_3.npy')

cd_4 = np.array(pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_4.pk1'))
sf_4 =  np.load('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_4.npy')

cd_5 = np.array(pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_5.pk1'))
sf_5 = np.load('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_5.npy')




combined_df =np.zeros(np.shape(cd_1), dtype = np.float32)
a,b = np.shape(cd_1)
for j in range(a):
    #Find in which band source is brightest 
    fluxArr = [-9999, cd_1[j,3], cd_2[j,3], cd_3[j,3], cd_4[j,3], cd_5[j,3] ]
    filt = np.argmax(fluxArr)
    #print (filt,cd_1[j,3], cd_2[j,3], cd_3[j,3], cd_4[j,3], cd_5[j,3])
    if(filt == 1):
        selected_coadd = cd_1
        selected_sf = sf_1
    elif(filt == 2):
        selected_coadd = cd_2
        selected_sf = sf_2
    elif(filt == 3):
        selected_coadd = cd_3
        selected_sf = sf_3
    elif(filt == 4):
        selected_coadd = cd_4
        selected_sf = sf_4
    elif(filt == 5):
        selected_coadd = cd_5
        selected_sf = sf_5
    
    wtArr = []
    sigxxArr=[]
    sigyyArr=[]
    sigxyArr=[]
    
    if(selected_coadd[j, 3]> 5):
        wtArr.append(200000)
        sigxxArr.append(selected_coadd[j, 7]  )
        sigyyArr.append(selected_coadd[j, 8] -selected_coadd[j, 39] )
        sigxyArr.append(selected_coadd[j, 9] -selected_coadd[j, 40])
    
# =============================================================================
#     for k in range(20):
#         #Measurement converged
#         if(selected_sf[k,j, 12] == 99 and selected_sf[k,j, 3]>0 and (not np.isnan(selected_sf[k,j, 7])) and selected_sf[k,j, 38] != -99):
#             wtArr.append(1)
#             #print (selected_sf[k,j, 7], selected_sf[k,j, 38])
#             sigxxArr.append(selected_sf[k,j, 7] )
#             sigyyArr.append(selected_sf[k,j, 8] )
#             sigxyArr.append(selected_sf[k,j, 9] )
#         #measurement forced
#         elif(selected_sf[k, j, 12] == 1 and selected_sf[k,j, 31]>0 and (not np.isnan(selected_sf[k,j, 35])) and selected_sf[k,j, 38]!= -99):
#             wtArr.append(0.5)
#             #print (selected_sf[k,j, 35], selected_sf[k,j, 38])
#             sigxxArr.append(selected_sf[k,j, 35] )
#             sigyyArr.append(selected_sf[k,j, 36] )
#             sigxyArr.append(selected_sf[k,j, 37] )
#         else:
#             wtArr.append(0.00001)
#             sigxxArr.append(1)
#             sigyyArr.append(1)
#             sigxyArr.append(1)
# =============================================================================
            
        combined_df[j, 7] = np.average(sigxxArr, weights= wtArr)
        combined_df[j, 8] = np.average(sigyyArr, weights= wtArr)
        combined_df[j, 9] = np.average(sigxyArr, weights= wtArr)
       
        
for j in range(mod_ySize):
    print (int(j*chopSize + chopSize/2))
    for k in range(mod_xSize):
        
        x_mid = int(k*chopSize + chopSize/2)
        y_mid = int(j*chopSize + chopSize/2)
        cond = np.where((ir_coadd_data[:,10] > x_mid-3000) & (ir_coadd_data[:,10] < x_mid+3000) & 
                        (ir_coadd_data[:,11] > y_mid-3000) & (ir_coadd_data[:,11] < y_mid+3000) 
                        & (ir_coadd_data[:,3] > 5) & (zArr > 0.2) & (ir_coadd_data[:,2] != 1)) 
        
        temp1 = np.copy(ir_coadd_data[cond])
        dx = temp1[:,10] -x_mid
        dy = temp1[:,11] -y_mid
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        
        #e1 = (temp1[:,7]-temp1[:,8])/(temp1[:,7]+temp1[:,8])
        #e2 = 2*temp1[:,9]/(temp1[:,7]+temp1[:,8])
        
        e1 = (combined_df[cond[0],7]-combined_df[cond[0],8])/(combined_df[cond[0],7]+combined_df[cond[0],8])
        e2 = 2*combined_df[cond[0],9]/(combined_df[cond[0],7]+combined_df[cond[0],8])
        
        
        epar= - (e1*cos2phi+e2*sin2phi)
        eper= (e2*cos2phi-e1*sin2phi)
        goodEllipIndices = np.where((np.abs(e1)<1) & (np.abs(e2)<1) & (r2<(3*alphax)**2) & (r2>100))
        wt = np.exp(-(r2/(2*alphax**2) ))
        
        mean,median,std = sigma_clipped_stats(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2)
        e1sum = np.sum(epar[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        e2sum = np.sum(eper[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        count = np.sum(wt[goodEllipIndices[0]]**2 *(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2) )
        
            
        if(count == 0):
            e1sum = 0
            e2sum = 0
            sys.exit()
            continue
        count = np.sqrt(0.5*count)
        e1sum = e1sum/count
        e2sum = e2sum/count
        #print (count, e1sum)
        
        
        imgE[j, k] = e1sum
        imgB[j ,k] = len(e1[goodEllipIndices[0]])
        del temp1
        
hdu = fits.PrimaryHDU(imgE)  
hdu.writeto('/scratch/halstead/d/dutta26/lsst/EMode2.fits', overwrite=True)

hdu = fits.PrimaryHDU(imgB)  
hdu.writeto('/scratch/halstead/d/dutta26/lsst/BMode2.fits', overwrite=True)
            

            
            
            
        
    
    
    



    
    
    
    