#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 08:38:58 2022

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
import wquantiles


badLocs = np.load('/home/dutta26/codes/stripeLoc.npy')

#Read Redshifts
f=open('/home/dutta26/zphot.out')
content = f.readlines()
f.close()
redShiftArr=[]
for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        redShiftArr.append(float((content[j].split())[1]))
redShiftArr = np.ones(54141)*9      
        
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1')
ir_coadd_data = np.array(ir_coadd_data)
redShiftArr = np.array(redShiftArr)        
z_min = 0.3
z_max = 10
bandList =['g', 'r', 'i']

#g_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_g.npy')
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_i.npy')
#z_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_z.npy')

g_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_g.pk1'))
r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r.pk1'))
i_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_i.pk1'))
#z_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1'))

bkgIndices = np.where( (redShiftArr>z_min )& (redShiftArr<z_max)& (ir_coadd_data[:,3] > 0.1) & (ir_coadd_data[:,3] < 100000)
                      & (ir_coadd_data[:,2] == 0))[0]

temp_coadd = np.zeros((40,40))

a1,b1 = np.shape(ir_coadd_data)
master_frame = np.zeros((a1,10), dtype = np.float32)
cnt1 = cnt2=cnt3= cnt4=cnt5=cnt6=cnt7=cnt8=cnt9=0



#Find total background 
totBkg = 0
totScale = 0
for value in r_sf_df[:,0,15]:
    if(value > 0):
        totBkg += value
for value in i_sf_df[:,0,15]:
    if(value > 0):
        totBkg += value

#Find eff scale factor
for value in r_sf_df[:,0,30]:
    if(value > 0):
        totScale += 1/value
for value in i_sf_df[:,0,30]:
    if(value > 0):
        totScale += 1/value


print ('aa')


for j in range(a1):
    
#for j in [21407]:    
    print (j)
    
    if(j not in bkgIndices):
        continue
    
    sigxx_arr=[]
    sigyy_arr=[]
    sigxy_arr=[]
    
    psfsigxx_arr=[]
    psfsigyy_arr=[]
    psfsigxy_arr=[]
    
    psfsigxx_uncertain_arr=[]
    psfsigyy_uncertain_arr=[]
    psfsigxy_uncertain_arr=[]
    
    wtArr= []
    area_arr=[]
    B_arr=[]
    N_arr=[]
    indexArr=[]
    coadd_measurement_flag = 0.0
    
    #If super faint use coadd measurements
    if(ir_coadd_data [j,3] < 1):
        
        sigxx_arr.append(ir_coadd_data [j,35])
        sigyy_arr.append(ir_coadd_data [j,36])
        sigxy_arr.append(ir_coadd_data [j,37])
        
        psfsigxx_arr.append(ir_coadd_data [j,38])
        psfsigyy_arr.append(ir_coadd_data [j,39])
        psfsigxy_arr.append(ir_coadd_data [j,40])
        
        psfsigxx_uncertain_arr.append(ir_coadd_data [j,41])
        psfsigyy_uncertain_arr.append(ir_coadd_data [j,42])
        psfsigxy_uncertain_arr.append(ir_coadd_data [j,43])
        
        area = 2*np.pi*np.sqrt(ir_coadd_data [j,7]* ir_coadd_data [j,8] - ir_coadd_data [j,9]**2)
        B = totBkg
        N = ir_coadd_data [j,3] *totScale
        error = np.sqrt(1/N + (4*area*B)/N**2 )
        wtArr.append(1/error)
        
        area_arr.append(area)
        N_arr.append(N)
        B_arr.append(B)
        indexArr.append(j)
        coadd_measurement_flag = 1
    
        
    else:
        
        #First do r band 
        a,b,c = np.shape(r_sf_df)
        for k in range(a):
            
            
            if(r_sf_df[k,j,12] <= 0):
                continue
            if(r_sf_df[k,j,13] == 1 or r_sf_df[k,j,14] == 1):
                continue
            if(r_sf_df[k,j,38] <= 0 or r_sf_df[k,j,39] <= 0 ):
                continue
            
            if(r_sf_df[k,j,12] == 99):
                if(r_sf_df[k,j,7] < 0 or r_sf_df[k,j,7]>70 or r_sf_df[k,j,8] < 0 or r_sf_df[k,j,8]>70):
                    continue
                if(r_sf_df[k,j,7] == None or np.isnan(r_sf_df[k,j,7]) or r_sf_df[k,j,8] == None or np.isnan(r_sf_df[k,j,8])):
                    continue
                if(r_sf_df[k,j,3]<= 0 or r_sf_df[k,j,3]== None or np.isnan(r_sf_df[k,j,3]) ):
                    continue
                    
                
            if(r_sf_df[k,j,12] == 1):
                if(r_sf_df[k,j,35] < 0 or r_sf_df[k,j,35]>70 or r_sf_df[k,j,36] < 0 or r_sf_df[k,j,36]>70):
                    continue
                if(r_sf_df[k,j,35] == None or np.isnan(r_sf_df[k,j,35]) or r_sf_df[k,j,36] == None or np.isnan(r_sf_df[k,j,36])):
                    continue
                if(r_sf_df[k,j,31]<= 0 or r_sf_df[k,j,31]== None or np.isnan(r_sf_df[k,j,31]) ):
                    continue
            
                
            if(r_sf_df[k,j,12] == 99):
                sigxx_arr.append(r_sf_df[k,j,7])
                sigyy_arr.append(r_sf_df[k,j,8])
                sigxy_arr.append(r_sf_df[k,j,9])
                area = 2*np.pi*np.sqrt(r_sf_df[k,j,7]* r_sf_df[k,j,8] - r_sf_df[k,j,9]**2)
                #N =r_sf_df[k,j,3]
                N = r_coadd_df [j,3] / r_sf_df[k,j,30]
                B = r_sf_df[k,j,6]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
                
            elif(r_sf_df[k,j,12] == 1):
                sigxx_arr.append(r_sf_df[k,j,35])
                sigyy_arr.append(r_sf_df[k,j,36])
                sigxy_arr.append(r_sf_df[k,j,37])
                temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + r_sf_df[k,j,38]
                temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + r_sf_df[k,j,39]
                temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + r_sf_df[k,j,40]
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                #N =r_sf_df[k,j,31]
                N = r_coadd_df [j,3] / r_sf_df[k,j,30]
                B = r_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            else:
                continue
                
            psfsigxx_arr.append(r_sf_df[k,j,38])
            psfsigyy_arr.append(r_sf_df[k,j,39])
            psfsigxy_arr.append(r_sf_df[k,j,40])
            
            psfsigxx_uncertain_arr.append(r_sf_df[k,j,41])
            psfsigyy_uncertain_arr.append(r_sf_df[k,j,42])
            psfsigxy_uncertain_arr.append(r_sf_df[k,j,43])
            
            
            error = np.sqrt(1/N + (4*area*B)/N**2 )
            if(badLocs[int(r_sf_df[k,j,47]), int(r_sf_df[k,j,48]), int(r_sf_df[k,j,49])] == 0):
                wtArr.append(0)
                #print ('aa')
            else:
                wtArr.append(1/error)
                
         
            
        #Now do I band
        a,b,c = np.shape(i_sf_df)
        for k in range(a):
            
            
            if(i_sf_df[k,j,12] <= 0):
                continue
            if(i_sf_df[k,j,13] == 1 or i_sf_df[k,j,14] == 1):
                continue
            if(i_sf_df[k,j,38] <= 0 or i_sf_df[k,j,39] <= 0 ):
                continue
            
            if(i_sf_df[k,j,12] == 99):
                if(i_sf_df[k,j,7] < 0 or i_sf_df[k,j,7]>70 or i_sf_df[k,j,8] < 0 or i_sf_df[k,j,8]>70):
                    continue
                if(i_sf_df[k,j,7] == None or np.isnan(i_sf_df[k,j,7]) or i_sf_df[k,j,8] == None or np.isnan(i_sf_df[k,j,8])):
                    continue
                if(i_sf_df[k,j,3]<= 0 or i_sf_df[k,j,3]== None or np.isnan(i_sf_df[k,j,3]) ):
                    continue
                    
                
            if(i_sf_df[k,j,12] == 1):
                if(i_sf_df[k,j,35] < 0 or i_sf_df[k,j,35]>70 or i_sf_df[k,j,36] < 0 or i_sf_df[k,j,36]>70):
                    continue
                if(i_sf_df[k,j,35] == None or np.isnan(i_sf_df[k,j,35]) or i_sf_df[k,j,36] == None or np.isnan(i_sf_df[k,j,36])):
                    continue
                if(i_sf_df[k,j,31]<= 0 or i_sf_df[k,j,31]== None or np.isnan(i_sf_df[k,j,31]) ):
                    continue
            
                
            if(i_sf_df[k,j,12] == 99):
                sigxx_arr.append(i_sf_df[k,j,7])
                sigyy_arr.append(i_sf_df[k,j,8])
                sigxy_arr.append(i_sf_df[k,j,9])
                area = 2*np.pi*np.sqrt(i_sf_df[k,j,7]* i_sf_df[k,j,8] - i_sf_df[k,j,9]**2)
                #N =i_sf_df[k,j,3]
                N = i_coadd_df [j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,6]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
                
            elif(i_sf_df[k,j,12] == 1):
                sigxx_arr.append(i_sf_df[k,j,35])
                sigyy_arr.append(i_sf_df[k,j,36])
                sigxy_arr.append(i_sf_df[k,j,37])
                temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + i_sf_df[k,j,38]
                temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + i_sf_df[k,j,39]
                temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + i_sf_df[k,j,40]
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                #N =i_sf_df[k,j,31]
                N = i_coadd_df[j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            else:
                continue
                
            psfsigxx_arr.append(i_sf_df[k,j,38])
            psfsigyy_arr.append(i_sf_df[k,j,39])
            psfsigxy_arr.append(i_sf_df[k,j,40])
            
            psfsigxx_uncertain_arr.append(i_sf_df[k,j,41])
            psfsigyy_uncertain_arr.append(i_sf_df[k,j,42])
            psfsigxy_uncertain_arr.append(i_sf_df[k,j,43])
            
            
            error = np.sqrt(1/N + (4*area*B)/N**2 )
            if(badLocs[int(i_sf_df[k,j,47]), int(i_sf_df[k,j,48]), int(i_sf_df[k,j,49])] == 0):
                wtArr.append(0)
                #print ('aa')
            else:
                wtArr.append(1/error)
                
        
        

                
    print (j,len(wtArr), len(indexArr))
    if(len(wtArr) == 0):
        continue
    #Now correct for PSF and use monte carlo if needed
    corr_sigxx_arr=[]
    corr_sigyy_arr=[]
    corr_sigxy_arr=[]
    for k in range(len(sigxx_arr)):
        corr_xx = sigxx_arr[k] - psfsigxx_arr[k]
        corr_yy = sigyy_arr[k] - psfsigyy_arr[k]
        corr_xy = sigxy_arr[k] - psfsigxy_arr[k]
        temp = corr_xx +corr_yy - 2*np.abs(corr_xy)
        
        if(corr_xx< 0 or corr_yy<0 or temp<0):
            
# =============================================================================
#             
#             
#             temp_xx = ir_coadd_data[j, 7] - ir_coadd_data[j, 38] + psfsigxx_arr[k]
#             temp_yy = ir_coadd_data[j, 8] - ir_coadd_data[j, 39] + psfsigyy_arr[k]
#             temp_xy = ir_coadd_data[j, 9] - ir_coadd_data[j, 40] + psfsigxy_arr[k]
#             e1_measured = (temp_xx - temp_yy)/(temp_xx +temp_yy)
#             e2_measured = 2*temp_xy/(temp_xx +temp_yy)
#                 
#             s = np.sqrt( (area_arr[k]/(np.pi*N_arr[k]) + 4*area_arr[k]**2 * B_arr[k]/(np.pi * N_arr[k]**2)) ) * area_arr[k]
#             s_xx = (1+e1_measured)*s
#             s_yy = (1-e1_measured)*s
#             s_xy = abs(e2_measured)*s
#             
#           
#             
#             corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx_arr[k], sigyy_arr[k], sigxy_arr[k], 
#                                                                  psfsigxx_arr[k], psfsigyy_arr[k], psfsigxy_arr[k], 
#                                                                  s_xx,s_yy,s_xy, 
#                                                                  psfsigxx_uncertain_arr[k], 
#                                                                  psfsigyy_uncertain_arr[k], 
#                                                                  psfsigxy_uncertain_arr[k])
#             
#             if(abs(corr_xx) > 70 or abs(corr_yy) >70 or abs(corr_xy)>70):
#                 corr_sigxx_arr.append(None)
#                 corr_sigyy_arr.append(None)
#                 corr_sigxy_arr.append(None)
#             elif(success < 50):
#                 corr_sigxx_arr.append(None)
#                 corr_sigyy_arr.append(None)
#                 corr_sigxy_arr.append(None)
#             else:
# =============================================================================
            corr_sigxx_arr.append(None)
            corr_sigyy_arr.append(None)
            corr_sigxy_arr.append(None)
        
        else:
            corr_sigxx_arr.append(corr_xx)
            corr_sigyy_arr.append(corr_yy)
            corr_sigxy_arr.append(corr_xy)
        
    e1_arr=[]
    e2_arr=[]
    for k in range(len(corr_sigxx_arr)):
        if(corr_sigxx_arr[k] == None  or np.isnan(corr_sigxx_arr[k]) or corr_sigxx_arr[k] <0):
            e1_arr.append(0)
            e2_arr.append(0)
            wtArr[k] = 0
            corr_sigxx_arr[k] = 0.0
            corr_sigyy_arr[k] = 0.0
            corr_sigxy_arr[k] = 0.0
            continue
        if(corr_sigyy_arr[k] == None or np.isnan(corr_sigyy_arr[k]) or corr_sigyy_arr[k] <0 ):
            e1_arr.append(0)
            e2_arr.append(0)
            wtArr[k] = 0
            corr_sigxx_arr[k] = 0.0
            corr_sigyy_arr[k] = 0.0
            corr_sigxy_arr[k] = 0.0
            continue
        
        if(corr_sigxy_arr[k] == None or np.isnan(corr_sigxy_arr[k])):
            e1_arr.append(0)
            e2_arr.append(0)
            wtArr[k] = 0
            corr_sigxx_arr[k] = 0.0
            corr_sigyy_arr[k] = 0.0
            corr_sigxy_arr[k] = 0.0
            continue
        
        e1 = (corr_sigxx_arr[k] - corr_sigyy_arr[k])/(corr_sigxx_arr[k] + corr_sigyy_arr[k])
        e2 = 2*corr_sigxy_arr[k] /(corr_sigxx_arr[k] + corr_sigyy_arr[k])
        
        e1_arr.append(e1)
        e2_arr.append(e2)
        
        
    e1_arr=np.array(e1_arr)
    e2_arr=np.array(e2_arr)
    wtArr = np.array(wtArr)
    
    #Combine optimally 
    e1 = np.sum(e1_arr*wtArr)/np.sum(wtArr)
    e2 = np.sum(e2_arr*wtArr)/np.sum(wtArr)
    if(np.isnan(e1) or np.isnan(e2)):
        sys.exit()
        e1 = e2 = 0.0
    #e1= wquantiles.median(e1_arr,wtArr)
    #e2= wquantiles.median(e2_arr,wtArr)
    
    
    
# =============================================================================
#     corr_sigxx_arr = np.array(corr_sigxx_arr)
#     corr_sigyy_arr = np.array(corr_sigyy_arr)
#     corr_sigxy_arr = np.array(corr_sigxy_arr)
#     wtArr = np.array(wtArr)
#     sigxx = np.sum(corr_sigxx_arr * wtArr)/ np.sum(wtArr)
#     sigyy = np.sum(corr_sigyy_arr * wtArr)/ np.sum(wtArr)
#     sigxy = np.sum(corr_sigxy_arr * wtArr)/ np.sum(wtArr)
#     e1 = (sigxx-sigyy)/(sigxx+sigyy)
#     e2 = 2*sigxy/(sigxx+sigyy)
# =============================================================================
    
    validIndices = np.where(wtArr >0)[0]
    
    master_frame[j, 0] = ir_coadd_data[j,10]
    master_frame[j, 1] = ir_coadd_data[j,11]
    master_frame[j, 2] = e1
    master_frame[j, 3] = e2
    master_frame[j, 4] = redShiftArr[j]
    master_frame[j,5] = len(validIndices)
    
    
# =============================================================================
#     if(len(wtArr) > 50 and j>5000 and ir_coadd_data[j,3]> 400 ):
#         print ('aa')
#         sys.exit()
# =============================================================================
#sys.exit()
np.save('/scratch/bell/dutta26/abell_2390/master_arr.npy', master_frame)
#sys.exit()




















master_frame = np.load('/scratch/bell/dutta26/abell_2390/master_arr.npy')
#Make the make wrt to ir coadd 
coadd_file = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits'
f=fits.open(coadd_file)
data = f[0].data
ySize,xSize = np.shape(data)
f.close()  
del data 
    
chopSize = 50
alphax = 1000

mod_ySize = int(round(ySize/chopSize)) + 1
mod_xSize = int(round(xSize/chopSize)) + 1
imgE = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
imgB = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
count_img = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)

for j in range(mod_ySize):
    print (int(j*chopSize + chopSize/2))
    for k in range(mod_xSize):
        
        x_mid = int(k*chopSize + chopSize/2)
        y_mid = int(j*chopSize + chopSize/2)
        cond = np.where((master_frame[:,0] > x_mid-3000) & (master_frame[:,0] < x_mid+3000) & 
                        (master_frame[:,1] > y_mid-3000) & (master_frame[:,1] < y_mid+3000) 
                        & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,2] !=0)  ) 
        

        temp = np.copy(master_frame[cond])
        dx = temp[:,0] -x_mid
        dy = temp[:,1] -y_mid
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        
        e1 = temp[:,2]
        e2 = temp[:,3]
        tot_ellip = np.sqrt(e1**2 + e2**2)
        epar= - (e1*cos2phi+e2*sin2phi)
        eper= (e2*cos2phi-e1*sin2phi)
        #goodEllipIndices = np.where((np.abs(e1)<0.8) & (np.abs(e2)<0.8) & (r2<(3*alphax)**2) & (r2>100))
        goodEllipIndices = np.where((tot_ellip < 0.8)  & (r2<(3*alphax)**2) & (r2>100))
        wt = np.exp(-(r2/(2*alphax**2) ))
        
        mean,median,std = sigma_clipped_stats(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2)
        e1sum = np.sum(epar[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        e2sum = np.sum(eper[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        count = np.sum(wt[goodEllipIndices[0]]**2 *(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2) )
        
        if(count == 0):
            e1sum = 0
            e2sum = 0
            #sys.exit()
            #continue
        count = np.sqrt(0.5*count)
        e1sum = e1sum/count
        e2sum = e2sum/count
        #print (count, e1sum, e2sum)
        
        
        imgE[j, k] = e1sum
        imgB[j ,k] = e2sum
        count_img[j ,k] = len(epar[goodEllipIndices[0]])
        del temp
        
hdu = fits.PrimaryHDU(imgE)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/EMode.fits', overwrite=True)

hdu = fits.PrimaryHDU(imgB)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/BMode.fits', overwrite=True)

hdu = fits.PrimaryHDU(count_img)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/count_img.fits', overwrite=True)
            

            
