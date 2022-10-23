#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:48:09 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import wquantiles,sys,os,shutil
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
from scipy.special import erf

ir_coadd_data = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_ir.pk1')
ir_coadd_data = np.array(ir_coadd_data)
zArr = np.load('/scratch/halstead/d/dutta26/lsst/lsst_zArr.npy')
min_dist_arr = np.load('/scratch/halstead/d/dutta26/lsst/min_dist_arr.npy')


ir_coadd_data[:,7] = ir_coadd_data[:,7] 
ir_coadd_data[:,8] = ir_coadd_data[:,8] 
ir_coadd_data[:,9] = ir_coadd_data[:,9] 
# =============================================================================
# for j in range(len(min_dist_arr)):
#     size =np.pi * np.sqrt(ir_coadd_data[j,7]+ir_coadd_data[j,8]) 
#     if((ir_coadd_data[j,7] -ir_coadd_data[j,38])< 0.15 ):
#         N = ir_coadd_data[j,3]
#         err_sigxx = np.sqrt((size/np.pi)/N + (4*size* 1500* size/(np.pi*N*N)))
#         
#         s = np.sqrt(err_sigxx**2 + 0.1**2 )
#         A = ir_coadd_data[j,7]
#         B = ir_coadd_data[j,38]
#         ir_coadd_data[j,7] = 0.5*(A-B)*(erf((A-B)/np.sqrt(2)/s) + 1)   + s/np.sqrt(2*np.pi)*np.exp(-(A-B)**2/2/s/s)
#     else:
#         ir_coadd_data[j,7]= (ir_coadd_data[j,7] -ir_coadd_data[j,38])
#         
#     if((ir_coadd_data[j,8] -ir_coadd_data[j,39])< 0.15 ):
#         N = ir_coadd_data[j,3]        
#         err_sigxx = np.sqrt((size/np.pi)/N + (4*size* 1500* size/(np.pi*N*N)))
#         
#         s = np.sqrt(err_sigxx**2 + 0.1**2 )
#         A = ir_coadd_data[j,8]
#         B = ir_coadd_data[j,39]
#         ir_coadd_data[j,8] = 0.5*(A-B)*(erf((A-B)/np.sqrt(2)/s) + 1)   + s/np.sqrt(2*np.pi)*np.exp(-(A-B)**2/2/s/s)
#     else:
#         ir_coadd_data[j,8]= (ir_coadd_data[j,8] -ir_coadd_data[j,39])
# =============================================================================
            

e1_temp = (ir_coadd_data[:,7]-ir_coadd_data[:,8])/(ir_coadd_data[:,7]+ir_coadd_data[:,8])
e2_temp = 2*ir_coadd_data[:,9]/(ir_coadd_data[:,7]+ir_coadd_data[:,8])
       
coadd_file = '/scratch/halstead/d/dutta26/lsst/coadd1_ir.fits'
f=fits.open(coadd_file)
data = np.array(f[0].data)
ySize,xSize = np.shape(data)
f.close()  
del data 
 

#Crete test1 and 2
os.remove('/scratch/halstead/d/dutta26/lsst/test1.fits')
os.remove('/scratch/halstead/d/dutta26/lsst/test.fits')
shutil.copy(coadd_file,'/scratch/halstead/d/dutta26/lsst/test1.fits')
shutil.copy(coadd_file,'/scratch/halstead/d/dutta26/lsst/test.fits')
# =============================================================================
# count = 0 
# ellipArr=[]
# f_test=fits.open('/scratch/halstead/d/dutta26/lsst/test.fits', mode = 'update')
# for j in range(len(zArr)):
#     if(ir_coadd_data[j,3] < 300 or zArr[j]< 0.2 or min_dist_arr[j]< 10 or ir_coadd_data[j,2] == 1 ):
#         continue
#     else:
#         x= int(ir_coadd_data[j,10])
#         y= int(ir_coadd_data[j,11])
#         if(x<25 or x>8900 or y<25 or y>8900):
#             continue
#         #if(abs(e1_temp[j])>1 or abs(e2_temp[j])>1):
#         #    ellipArr.append(j)
#         if(ir_coadd_data[j,7]<0 or ir_coadd_data[j,8]< 0):
#             ellipArr.append(np.sqrt(e1_temp[j]**2 + e2_temp[j]**2))
#             continue
#         f_test[0].data[y-20, x-20:x+20] = 1000
#         f_test[0].data[y+20, x-20:x+20] = 1000
#         f_test[0].data[y-20:y+20, x-20] = 1000
#         f_test[0].data[y-20:y+20, x+20] = 1000 
#         count += 1
# 
# 
# f_test.flush()
# sys.exit()
# =============================================================================
   
chopSize = 25
alphax = 250

mod_ySize = int(round(ySize/chopSize)) + 2
mod_xSize = int(round(xSize/chopSize)) + 2
imgE = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
imgB = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)


#master_i = master_frame[0,:,:]
#master_r = master_frame[1,:,:]
        
        
for j in range(mod_ySize):
    print (int(j*chopSize + chopSize/2))
    for k in range(mod_xSize):
        
        x_mid = int(k*chopSize + chopSize/2)
        y_mid = int(j*chopSize + chopSize/2)
        cond = np.where((ir_coadd_data[:,10] > x_mid-1000) & (ir_coadd_data[:,10] < x_mid+1000) & 
                        (ir_coadd_data[:,11] > y_mid-1000) & (ir_coadd_data[:,11] < y_mid+1000) 
                        & (ir_coadd_data[:,3] > 300) & (zArr > 0.2) & (min_dist_arr > 10)) 
        
        temp1 = np.copy(ir_coadd_data[cond])
        dx = temp1[:,10] -x_mid
        dy = temp1[:,11] -y_mid
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        
        e1 = (temp1[:,7]-temp1[:,8])/(temp1[:,7]+temp1[:,8])
        e2 = 2*temp1[:,9]/(temp1[:,7]+temp1[:,8])
        
        epar= - (e1*cos2phi+e2*sin2phi)
        eper= (e2*cos2phi-e1*sin2phi)
        goodEllipIndices = np.where((np.abs(e1)<0.99) & (np.abs(e2)<0.99) & (r2<(3*alphax)**2) & (r2>100))
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
        imgB[j ,k] = e2sum
        del temp1
        
hdu = fits.PrimaryHDU(imgE)  
hdu.writeto('/scratch/halstead/d/dutta26/lsst/EMode2.fits', overwrite=True)

hdu = fits.PrimaryHDU(imgB)  
hdu.writeto('/scratch/halstead/d/dutta26/lsst/BMode2.fits', overwrite=True)
            

            
            
            
print (np.sum(imgE[175:185, 175:185])  )      
    
#504.66602    
    
#609 current record with no correction and only birght(300) seprated(10) objects 
#604 is reocrd using bright stars and correction + 1.15

    
    
    
    