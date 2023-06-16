#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 11:24:45 2022

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


f=open('/home/dutta26/zphot.out')
content = f.readlines()
f.close()

redShiftArr=[]

for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        redShiftArr.append(float((content[j].split())[1]))
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
ir_coadd_data = np.array(ir_coadd_data)
redShiftArr = np.array(redShiftArr)        
z_min = 0
z_max = 10
bandList =['g', 'r', 'i', 'z']

#g_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_g.npy')
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test4t3_r.npy')
i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test4t3_i.npy')
#z_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_z.npy')

#g_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_g_.pk1'))
r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r_.pk1'))
i_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_i_.pk1'))
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

e1_coadd_arr=[]
e1_sf_arr =[]

e2_coadd_arr =[]
e2_sf_arr=[]
j_arr =[]

measure_times =[]

for j in range(a1):
    

#for j in [32391]:
    
    if(j not in bkgIndices):
        continue
    
    if(ir_coadd_data [j,3] <0.75 or ir_coadd_data [j,3] >1.25):
        continue
    cnt2 += 1
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
    if(True):
        if(ir_coadd_data [j,35] == 0 ):
            sigxx_arr.append(ir_coadd_data [j,7])
            sigyy_arr.append(ir_coadd_data [j,8])
            sigxy_arr.append(ir_coadd_data [j,9])
        else:
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
    
        
    if(True):
        
        #First do r band 
        a,b,c = np.shape(r_sf_df)
        for k in range(a):
            
            
            if(r_sf_df[k,j,12] <= 0):
                continue
            if(r_sf_df[k,j,13] == 1 or r_sf_df[k,j,14] == 1 or r_sf_df[k,j,45] == 1 or r_sf_df[k,j,46] == 1):
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
                    
                
            elif(r_sf_df[k,j,12] == 1 and r_sf_df[k,j,50] == 0):
                if(r_sf_df[k,j,35] < 0 or r_sf_df[k,j,35]>70 or r_sf_df[k,j,36] < 0 or r_sf_df[k,j,36]>70):
                    continue
                if(r_sf_df[k,j,35] == None or np.isnan(r_sf_df[k,j,35]) or r_sf_df[k,j,36] == None or np.isnan(r_sf_df[k,j,36])):
                    continue
                if(r_sf_df[k,j,31]<= 0 or r_sf_df[k,j,31]== None or np.isnan(r_sf_df[k,j,31]) ):
                    continue
                
            elif(r_sf_df[k,j,12] == 1 and r_sf_df[k,j,50] != 0 and (not np.isnan(r_sf_df[k,j,50]))):
                if(r_sf_df[k,j,50] < 0 or r_sf_df[k,j,50]>70 or r_sf_df[k,j,51] < 0 or r_sf_df[k,j,51]>70):
                    continue
                if(r_sf_df[k,j,50] == None or np.isnan(r_sf_df[k,j,50]) or r_sf_df[k,j,51] == None or np.isnan(r_sf_df[k,j,51])):
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
                
            elif(r_sf_df[k,j,12] == 1 and (r_sf_df[k,j,50] == 0 or np.isnan(r_sf_df[k,j,50]))):
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
            elif(r_sf_df[k,j,12] == 1 and r_sf_df[k,j,50] != 0 and (not np.isnan(r_sf_df[k,j,50]))):
                sigxx_arr.append(r_sf_df[k,j,50] +r_sf_df[k,j,38])
                sigyy_arr.append(r_sf_df[k,j,51]+ r_sf_df[k,j,39])
                sigxy_arr.append(r_sf_df[k,j,52] + r_sf_df[k,j,40] )
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
                #print ('aa')
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
            if(i_sf_df[k,j,13] == 1 or i_sf_df[k,j,14] == 1 or i_sf_df[k,j,45] == 1 or i_sf_df[k,j,46] == 1):
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
                    
                
            if(i_sf_df[k,j,12] == 1 and i_sf_df[k,j,50] == 0):
                if(i_sf_df[k,j,35] < 0 or i_sf_df[k,j,35]>70 or i_sf_df[k,j,36] < 0 or i_sf_df[k,j,36]>70):
                    continue
                if(i_sf_df[k,j,35] == None or np.isnan(i_sf_df[k,j,35]) or i_sf_df[k,j,36] == None or np.isnan(i_sf_df[k,j,36])):
                    continue
                if(i_sf_df[k,j,31]<= 0 or i_sf_df[k,j,31]== None or np.isnan(i_sf_df[k,j,31]) ):
                    continue
                
            elif(i_sf_df[k,j,12] == 1 and i_sf_df[k,j,50] != 0 and (not np.isnan(i_sf_df[k,j,50]))):
                if(i_sf_df[k,j,50] < 0 or i_sf_df[k,j,50]>70 or i_sf_df[k,j,51] < 0 or i_sf_df[k,j,51]>70):
                    continue
                if(i_sf_df[k,j,50] == None or np.isnan(i_sf_df[k,j,50]) or i_sf_df[k,j,51] == None or np.isnan(i_sf_df[k,j,51])):
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
                
            elif(i_sf_df[k,j,12] == 1 and (i_sf_df[k,j,50] == 0 or np.isnan(i_sf_df[k,j,50]))):
                sigxx_arr.append(i_sf_df[k,j,35])
                sigyy_arr.append(i_sf_df[k,j,36])
                sigxy_arr.append(i_sf_df[k,j,37])
                temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + i_sf_df[k,j,38]
                temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + i_sf_df[k,j,39]
                temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + i_sf_df[k,j,40]
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                #N =i_sf_df[k,j,31]
                N = i_coadd_df [j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            elif(i_sf_df[k,j,12] == 1 and i_sf_df[k,j,50] != 0 and (not np.isnan(i_sf_df[k,j,50]))):
                sigxx_arr.append(i_sf_df[k,j,50] +i_sf_df[k,j,38])
                sigyy_arr.append(i_sf_df[k,j,51]+ i_sf_df[k,j,39])
                sigxy_arr.append(i_sf_df[k,j,52] + i_sf_df[k,j,40] )
                temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + i_sf_df[k,j,38]
                temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + i_sf_df[k,j,39]
                temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + i_sf_df[k,j,40]
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                #N =i_sf_df[k,j,31]
                N = i_coadd_df [j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
                #print ('aa')
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
                
        
        

    #print (cnt1, cnt2, cnt3, cnt4, cnt5)            
    print (j,len(wtArr), len(indexArr), '**')
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
            
            
            
            temp_xx = ir_coadd_data[j, 7] - ir_coadd_data[j, 38] + psfsigxx_arr[k]
            temp_yy = ir_coadd_data[j, 8] - ir_coadd_data[j, 39] + psfsigyy_arr[k]
            temp_xy = ir_coadd_data[j, 9] - ir_coadd_data[j, 40] + psfsigxy_arr[k]
            e1_measured = (temp_xx - temp_yy)/(temp_xx +temp_yy)
            e2_measured = 2*temp_xy/(temp_xx +temp_yy)
                
            s = np.sqrt( (area_arr[k]/(np.pi*N_arr[k]) + 4*area_arr[k]**2 * B_arr[k]/(np.pi * N_arr[k]**2)) ) * area_arr[k]
            s_xx = (1+e1_measured)*s
            s_yy = (1-e1_measured)*s
            s_xy = abs(e2_measured)*s
            
          
            
            corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx_arr[k], sigyy_arr[k], sigxy_arr[k], 
                                                                 psfsigxx_arr[k], psfsigyy_arr[k], psfsigxy_arr[k], 
                                                                 s_xx,s_yy,s_xy, 
                                                                 psfsigxx_uncertain_arr[k], 
                                                                 psfsigyy_uncertain_arr[k], 
                                                                 psfsigxy_uncertain_arr[k])
            
            if(abs(corr_xx) > 70 or abs(corr_yy) >70 or abs(corr_xy)>70):
                corr_sigxx_arr.append(None)
                corr_sigyy_arr.append(None)
                corr_sigxy_arr.append(None)
            elif(success < 50):
                corr_sigxx_arr.append(None)
                corr_sigyy_arr.append(None)
                corr_sigxy_arr.append(None)
            else:
                corr_sigxx_arr.append(corr_xx)
                corr_sigyy_arr.append(corr_yy)
                corr_sigxy_arr.append(corr_xy)
        
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
    e1 = np.sum(e1_arr[1:]*wtArr[1:])/np.sum(wtArr[1:])
    e2 = np.sum(e2_arr[1:]*wtArr[1:])/np.sum(wtArr[1:])
    
    #e1= wquantiles.median(e1_arr,wtArr)
    #e2= wquantiles.median(e2_arr,wtArr)
    
    if(e1_arr[0] == 0 or e2_arr[0] ==0):
        #cnt1 += 1
        #j_arr.append(j)
        continue
    #if(abs(e1_arr[0]-e1)>0.5):
        #j_arr.append(j)
        #continue
    e1_coadd_arr.append(e1_arr[0])
    e1_sf_arr.append(e1)
    
    e2_coadd_arr.append(e2_arr[0])
    e2_sf_arr.append(e2)
    
    
    loc = np.where(wtArr>0)[0]
    measure_times.append(len(loc))


import matplotlib.pyplot as plt 
snr_coadd = 1*totScale / np.sqrt(1*totScale + 4*3.14*3*3*totBkg)
fig, axs = plt.subplots(2)
fig.suptitle('WithOUT Sub-Coadd. Coadd SNR~'+str(snr_coadd)[0:4])
axs[0].plot(e1_coadd_arr, e1_sf_arr, 'r.', markersize = 3)
axs[0].set(xlabel = 'e1_coadd', ylabel = 'e1_single_frames_combined')
axs[1].plot(e2_coadd_arr, e2_sf_arr, 'r.', markersize = 3)
axs[1].set(xlabel = 'e2_coadd', ylabel = 'e2_single_frames_combined')

# =============================================================================
# e1_coadd_arr = np.array(e1_coadd_arr)
# e2_coadd_arr = np.array(e2_coadd_arr)
# e1_sf_arr = np.array(e1_sf_arr)
# e2_sf_arr = np.array(e2_sf_arr)
# j_arr = np.array(j_arr)
# 
# loc = np.where((np.abs(e2_coadd_arr - e2_sf_arr) > 0.15) & (np.abs(e2_sf_arr ) < 0.1 ) 
#                & (np.abs(e2_coadd_arr ) > 0.1 ))[0]
# 
# 
# 
# 
# 
# 
# from shutil import copy
# copy('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits' , '/scratch/bell/dutta26/abell_2390/temp1.fits')
# 
# f=fits.open('/scratch/bell/dutta26/abell_2390/temp1.fits', mode='update')
# for indices in loc:
#     sf_index = j_arr[indices]
#     xLoc = int(ir_coadd_data[sf_index, 10])
#     yLoc = int(ir_coadd_data[sf_index, 11])
#     
#     f[0].data[yLoc-20:yLoc+20, xLoc-20] = 1000*sf_index
#     f[0].data[yLoc-20:yLoc+20, xLoc+20] = 1000*sf_index
#     f[0].data[yLoc-20,xLoc-20: xLoc+20] = 1000*sf_index
#     f[0].data[yLoc+20,xLoc-20: xLoc+20] = 1000*sf_index 
#     
# f.flush()
# 
# 
# =============================================================================


