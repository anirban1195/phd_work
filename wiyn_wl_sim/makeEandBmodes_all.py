#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 06:25:12 2023

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
from astropy import wcs
import matplotlib.pyplot as plt
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')


#f1 = fits.open('/scratch/bell/dutta26/abell_2390/temp.fits', mode = 'update')

def weighted_quantiles_interpolate(values, weights, quantiles=0.5):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, np.array(quantiles) * c[-1])]]


def EandB_all(ir_coadd_data_name, r_coadd_npy_name, i_coadd_npy_name, zFile, outFile, coadd_img, r_sf_name, i_sf_name):
    
    #redShiftArr = np.ones(15106)*9      
    r_sf_npy = np.load(r_sf_name)
    i_sf_npy = np.load(i_sf_name)
    r_coadd_npy = np.load(r_coadd_npy_name)
    i_coadd_npy = np.load(i_coadd_npy_name)
    #ir_coadd_data = pd.read_pickle(ir_coadd_data)
    #ir_coadd_data = np.array(ir_coadd_data)
    ir_coadd_data = np.load(ir_coadd_data_name)
    
    z_min = 0.4
    z_max = 10
    bandList =['g', 'r', 'i']
    if (zFile == None):
        redShiftArr = np.ones(len(ir_coadd_data))*9      
    else:
        #Read Redshifts
        f=open(zFile)
        content = f.readlines()
        f.close()
        redShiftArr=[]
        for j in range(len(content)):
            if (content[j][0] == '#'):
                continue
            else:
                redShiftArr.append(float((content[j].split())[1]))
    redShiftArr = np.array(redShiftArr)        
    #g_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_g.npy')
    #r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
    #i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_i.npy')
    #z_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_z.npy')
    
    #g_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_g.pk1'))
    #r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r.pk1'))
    #i_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_i.pk1'))
    #z_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1'))
    
    bkgIndices = np.where( (redShiftArr>z_min )& (redShiftArr<z_max)& (ir_coadd_data[:,3] > 10) & (ir_coadd_data[:,3] < 1e10)
                          & (ir_coadd_data[:,2] == 0))[0]
    print (len(bkgIndices))
    #sys.exit()
    temp_coadd = np.zeros((40,40))
    
    a1,b1 = np.shape(ir_coadd_data)
    master_frame = np.zeros((a1,20, 60), dtype = np.float32)
    cnt1 = cnt2=cnt3= cnt4=cnt5=cnt6=cnt7=cnt8=cnt9=0
    
    
    
    #Find total background 
    totBkg = 18000
    totScale = 800
    # =============================================================================
    # for value in r_sf_df[:,0,15]:
    #     if(value > 0):
    #         totBkg += value
    # for value in i_sf_df[:,0,15]:
    #     if(value > 0):
    #         totBkg += value
    # 
    # #Find eff scale factor
    # for value in r_sf_df[:,0,30]:
    #     if(value > 0):
    #         totScale += 1/value
    # for value in i_sf_df[:,0,30]:
    #     if(value > 0):
    #         totScale += 1/value
    # =============================================================================
    
    failedArr=[]
    print ('aa')
    arr1=[]
    arr2=[]
    arr3=[]
    arr_sf_wt =[]
    arr_sf_wt1 =[]
    for j in range(a1):
        
    #for j in [17420]:   
        #print (j)
        
        if(j not in bkgIndices or ir_coadd_data [j,3] == 0  or ir_coadd_data [j,13] == 1 or ir_coadd_data [j,14] == 1 or  ir_coadd_data [j,59] == 1):
            continue
        
        sigxx_arr=[]
        sigyy_arr=[]
        sigxy_arr=[]
        
        sigxx_err_arr=[]
        sigyy_err_arr=[]
        sigxy_err_arr=[]
        
        psfsigxx_arr=[]
        psfsigyy_arr=[]
        psfsigxy_arr=[]
        
        psfsigxx_uncertain_arr=[]
        psfsigyy_uncertain_arr=[]
        psfsigxy_uncertain_arr=[]
        
        size_sf_arr =[]
        
        wtArr= []
        area_arr=[]
        B_arr=[]
        N_arr=[]
        N_corr_arr=[]
        indexArr=[]
        coadd_measurement_flag = 0.0
    
        #If super faint use coadd measurements
        if(ir_coadd_data [j,3] < 10 or r_coadd_npy[j,3]<=0 or i_coadd_npy[j,3]<=0 ):
        
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
            size = np.sqrt(ir_coadd_data [j,7]+ ir_coadd_data [j,8])
            B = totBkg
            N = ir_coadd_data [j,3] *totScale
            e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
            error = ir_coadd_data[j,68] #np.sqrt((size**4/N + (4*np.pi*size**6 *B)/N**2)) 
            if(error <= 0 or error== None or np.isnan(error)):
                error = 100000000
                
            
            wtArr.append(1/error)
            
            area_arr.append(area)
            N_arr.append(N)
            B_arr.append(B)
            indexArr.append(j)
            coadd_measurement_flag = 1
            
            sigxx_err_arr.append( np.sqrt(ir_coadd_data [j,71]**2 + ir_coadd_data [j,41]**2 ))
            sigyy_err_arr.append( np.sqrt(ir_coadd_data [j,72]**2 + ir_coadd_data [j,42]**2 ))
            sigxy_err_arr.append( np.sqrt(ir_coadd_data [j,73]**2 + ir_coadd_data [j,43]**2 ))
            
            size_sf_arr.append(np.sqrt(ir_coadd_data [j,38] + ir_coadd_data [j,39]))
    
        else:
            
            #First do r band 
            a,b,c = np.shape(r_sf_npy)
            for k in range(a):
                
                if( abs((ir_coadd_data[j,7] - r_coadd_npy[j,7])/ir_coadd_data[j,7] ) >2  ):
                    continue
                if( abs((ir_coadd_data[j,8] - r_coadd_npy[j,8])/ir_coadd_data[j,8] ) >2  ):
                    continue
                if( r_coadd_npy[j,7]<=0 or r_coadd_npy[j,8]<=0):
                    continue
                
                if(r_sf_npy[k,j,12] <= 0 or r_sf_npy[k,j,13]> 0 or r_sf_npy[k,j,14]> 0 or r_sf_npy[k,j,38]==-99):
                    continue
                if(np.sum(r_sf_npy[k,j,61:67])>=1 or r_sf_npy[k,j,38]==-99):
                  
                    continue
                if(r_sf_npy[k,j,12] == 99):
                    if(r_sf_npy[k,j,7] < 0 or r_sf_npy[k,j,7]>70 or r_sf_npy[k,j,8] < 0 or r_sf_npy[k,j,8]>70):
                        continue
                    if(r_sf_npy[k,j,7] == None or np.isnan(r_sf_npy[k,j,7]) or r_sf_npy[k,j,8] == None or np.isnan(r_sf_npy[k,j,8])):
                        continue
                    if(r_sf_npy[k,j,3]<= 0 or r_sf_npy[k,j,3]== None or np.isnan(r_sf_npy[k,j,3]) ):
                        continue
                        
                    
                if(r_sf_npy[k,j,12] == 1):
                    if(r_sf_npy[k,j,35] < 0 or r_sf_npy[k,j,35]>70 or r_sf_npy[k,j,36] < 0 or r_sf_npy[k,j,36]>70):
                        continue
                    if(r_sf_npy[k,j,35] == None or np.isnan(r_sf_npy[k,j,35]) or r_sf_npy[k,j,36] == None or np.isnan(r_sf_npy[k,j,36])):
                        continue
                    if(r_sf_npy[k,j,31]<= 0 or r_sf_npy[k,j,31]== None or np.isnan(r_sf_npy[k,j,31]) ):
                        
                        continue
                
                
                if(r_sf_npy[k,j,12] == 99):
                    sigxx_arr.append(r_sf_npy[k,j,51])
                    sigyy_arr.append(r_sf_npy[k,j,52])
                    sigxy_arr.append(r_sf_npy[k,j,53])
                    area = 2*np.pi*np.sqrt(r_sf_npy[k,j,7]* r_sf_npy[k,j,8] - r_sf_npy[k,j,9]**2)
                    size = np.sqrt(r_sf_npy[k,j,7]+r_sf_npy[k,j,8])
                    N =r_sf_npy[k,j,3]
                    #N = r_coadd_df [j,3] / r_sf_npy[k,j,30]
                    B = r_sf_npy[k,j,6]
                    area_arr.append(area)
                    if(N< 0 or N>1e5 ):
                        N = 0.00001
                    N_arr.append(N)
                    N_corr_arr.append(N*r_sf_npy[k,j,30])
                    B_arr.append(B)
                    indexArr.append(j)
                   
                    
                elif(r_sf_npy[k,j,12] == 1):
                    sigxx_arr.append(r_sf_npy[k,j,51])
                    sigyy_arr.append(r_sf_npy[k,j,52])
                    sigxy_arr.append(r_sf_npy[k,j,53])
                    temp_xx = r_sf_npy[k,j,51]
                    temp_yy = r_sf_npy[k,j,52]
                    temp_xy = r_sf_npy[k,j,53]
                    area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                    size = np.sqrt(temp_xx+temp_yy)
                    N =r_sf_npy[k,j,31]
                    if(N< 0 or N>1e5):
                        N = 0.00001
                    B = r_sf_npy[k,j,34]
                    area_arr.append(area)
                    N_arr.append(N)
                    N_corr_arr.append(N*r_sf_npy[k,j,30])
                    B_arr.append(B)
                    indexArr.append(j)
                    
                
                psfsigxx_arr.append(r_sf_npy[k,j,38])
                psfsigyy_arr.append(r_sf_npy[k,j,39])
                psfsigxy_arr.append(r_sf_npy[k,j,40])
                
                psfsigxx_uncertain_arr.append(r_sf_npy[k,j,41])
                psfsigyy_uncertain_arr.append(r_sf_npy[k,j,42])
                psfsigxy_uncertain_arr.append(r_sf_npy[k,j,43])
                
                sigxx_err_arr.append( np.sqrt(r_sf_npy[k,j,71]**2 + r_sf_npy[k,j,41]**2 ))
                sigyy_err_arr.append( np.sqrt(r_sf_npy[k,j,72]**2 + r_sf_npy[k,j,42]**2 ))
                sigxy_err_arr.append( np.sqrt(r_sf_npy[k,j,73]**2 + r_sf_npy[k,j,43]**2 ))
                size_sf_arr.append(np.sqrt(r_sf_npy[k,j,38] + r_sf_npy[k,j,39]))
                
                
                e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
                e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
                e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
                
                turb_xx = np.sqrt(r_sf_npy[k,j,41]**2 )
                turb_yy = np.sqrt(r_sf_npy[k,j,42]**2 )
                turb_avg = 0.5*(turb_xx+turb_yy)
                
               
                
                
                error =np.sqrt(size**4/N + (4*np.pi*size**6*B)/N**2 +turb_avg)
                if(error == 0 or error == None or np.isnan(error)):
                    error = 10000000000
                    wtArr.append(1/error)
                else:
                    
                    wtArr.append(1/error)
                    
             
            print (len(wtArr))   
            #Now do I band
            a,b,c = np.shape(i_sf_npy)
            for k in range(a):
                
                if( abs((ir_coadd_data[j,7] - i_coadd_npy[j,7])/ir_coadd_data[j,7] ) >2  ):
                    continue
                if( abs((ir_coadd_data[j,8] - i_coadd_npy[j,8])/ir_coadd_data[j,8] ) >2  ):
                    continue
                if( i_coadd_npy[j,7]<=0 or i_coadd_npy[j,8]<=0):
                    continue
                if(i_sf_npy[k,j,12] <= 0 or i_sf_npy[k,j,13]> 0 or i_sf_npy[k,j,14]> 0) or i_sf_npy[k,j,38]==-99:
                    continue
                if(np.sum(i_sf_npy[k,j,61:67])>=1 or i_sf_npy[k,j,38]==-99):
                    continue
                if(i_sf_npy[k,j,12] == 99):
                    if(i_sf_npy[k,j,7] < 0 or i_sf_npy[k,j,7]>70 or i_sf_npy[k,j,8] < 0 or i_sf_npy[k,j,8]>70):
                        continue
                    if(i_sf_npy[k,j,7] == None or np.isnan(i_sf_npy[k,j,7]) or i_sf_npy[k,j,8] == None or np.isnan(i_sf_npy[k,j,8])):
                        continue
                    if(i_sf_npy[k,j,3]<= 0 or i_sf_npy[k,j,3]== None or np.isnan(i_sf_npy[k,j,3]) ):
                        continue
                        
                    
                if(i_sf_npy[k,j,12] == 1):
                    if(i_sf_npy[k,j,35] < 0 or i_sf_npy[k,j,35]>70 or i_sf_npy[k,j,36] < 0 or i_sf_npy[k,j,36]>70):
                        continue
                    if(i_sf_npy[k,j,35] == None or np.isnan(i_sf_npy[k,j,35]) or i_sf_npy[k,j,36] == None or np.isnan(i_sf_npy[k,j,36])):
                        continue
                    if(i_sf_npy[k,j,31]<= 0 or i_sf_npy[k,j,31]== None or np.isnan(i_sf_npy[k,j,31]) ):
                        continue
                
                
                if(i_sf_npy[k,j,12] == 99):
                    sigxx_arr.append(i_sf_npy[k,j,51])
                    sigyy_arr.append(i_sf_npy[k,j,52])
                    sigxy_arr.append(i_sf_npy[k,j,53])
                    area = 2*np.pi*np.sqrt(i_sf_npy[k,j,7]* i_sf_npy[k,j,8] - i_sf_npy[k,j,9]**2)
                    size = np.sqrt(i_sf_npy[k,j,7]+ i_sf_npy[k,j,8])
                    N =i_sf_npy[k,j,3]
                    #N = i_coadd_df [j,3] / i_sf_npy[k,j,30]
                    B = i_sf_npy[k,j,6]
                    area_arr.append(area)
                    if(N< 0 or N>1e5 ):
                        N = 0.00001
                    N_arr.append(N)
                    N_corr_arr.append(N*i_sf_npy[k,j,30])
                    B_arr.append(B)
                    indexArr.append(j)
                    
                   
                    
                elif(i_sf_npy[k,j,12] == 1):
                    sigxx_arr.append(i_sf_npy[k,j,51])
                    sigyy_arr.append(i_sf_npy[k,j,52])
                    sigxy_arr.append(i_sf_npy[k,j,53])
                    temp_xx = i_sf_npy[k,j,51] 
                    temp_yy = i_sf_npy[k,j,52] 
                    temp_xy = i_sf_npy[k,j,53] 
                    area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                    size = np.sqrt(temp_xx+temp_yy)
                    N = i_sf_npy[k,j,31]
                    if(N< 0 or N>1e5 ):
                        N = 0.00001
                    #N = i_coadd_df[j,3] / i_sf_npy[k,j,30]
                    B = i_sf_npy[k,j,34]
                    area_arr.append(area)
                    N_arr.append(N)
                    N_corr_arr.append(N*i_sf_npy[k,j,30])
                    B_arr.append(B)
                    indexArr.append(j)
                    
                    
                
                psfsigxx_arr.append(i_sf_npy[k,j,38])
                psfsigyy_arr.append(i_sf_npy[k,j,39])
                psfsigxy_arr.append(i_sf_npy[k,j,40])
                
                psfsigxx_uncertain_arr.append(i_sf_npy[k,j,41])
                psfsigyy_uncertain_arr.append(i_sf_npy[k,j,42])
                psfsigxy_uncertain_arr.append(i_sf_npy[k,j,43])
                
                sigxx_err_arr.append( np.sqrt(i_sf_npy[k,j,71]**2 + i_sf_npy[k,j,41]**2 ))
                sigyy_err_arr.append( np.sqrt(i_sf_npy[k,j,72]**2 + i_sf_npy[k,j,42]**2 ))
                sigxy_err_arr.append( np.sqrt(i_sf_npy[k,j,73]**2 + i_sf_npy[k,j,43]**2 ))
                size_sf_arr.append(np.sqrt(i_sf_npy[k,j,38] + i_sf_npy[k,j,39]))
                
                e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
                e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
                e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
                
                turb_xx = np.sqrt(i_sf_npy[k,j,41]**2 )
                turb_yy = np.sqrt(i_sf_npy[k,j,42]**2 )
                turb_avg = 0.5*(turb_xx+turb_yy)
                
                
                
                error = np.sqrt(size**4/N + (4*np.pi*size**6*B)/N**2 +turb_avg)
                
                if(error == 0 or error == None or np.isnan(error)):
                    error = 10000000000
                    wtArr.append(1/error)
                else:
                    wtArr.append(1/error)
                    
            
            
    
                    
        print (j,len(wtArr), len(indexArr))
        
        if(len(wtArr) <= 10):
            
            sigxx_arr=[]
            sigyy_arr=[]
            sigxy_arr=[]
            
            sigxx_err_arr=[]
            sigyy_err_arr=[]
            sigxy_err_arr=[]
            
            psfsigxx_arr=[]
            psfsigyy_arr=[]
            psfsigxy_arr=[]
            
            psfsigxx_uncertain_arr=[]
            psfsigyy_uncertain_arr=[]
            psfsigxy_uncertain_arr=[]
            
            size_sf_arr =[]
            
            wtArr= []
            area_arr=[]
            B_arr=[]
            N_arr=[]
            N_corr_arr=[]
            indexArr=[]
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
            size = np.sqrt(ir_coadd_data [j,7]+ ir_coadd_data [j,8])
            B = totBkg
            N = ir_coadd_data [j,3] *totScale
            e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
            turb_xx = np.sqrt(ir_coadd_data[j,41]**2 )
            turb_yy = np.sqrt(ir_coadd_data[j,42]**2 )
            turb_avg = 0.5*(turb_xx+turb_yy)
            
            
            error = np.sqrt(ir_coadd_data[j,68] + turb_avg) #np.sqrt((size**4/N + (4*np.pi*size**6 *B)/N**2)) 
            if(error <= 0 or error== None or np.isnan(error)):
                error = 100000000
                
            
            wtArr.append(1/error)
            
            area_arr.append(area)
            N_arr.append(N)
            B_arr.append(B)
            indexArr.append(j)
            coadd_measurement_flag = 1
            
            sigxx_err_arr.append( np.sqrt(ir_coadd_data [j,71]**2 + ir_coadd_data [j,41]**2 ))
            sigyy_err_arr.append( np.sqrt(ir_coadd_data [j,72]**2 + ir_coadd_data [j,42]**2 ))
            sigxy_err_arr.append( np.sqrt(ir_coadd_data [j,73]**2 + ir_coadd_data [j,43]**2 ))
            
            size_sf_arr.append(np.sqrt(ir_coadd_data [j,38] + ir_coadd_data [j,39]))
            
        if(len(wtArr)<=0):
            continue
        
        if(len(wtArr)== 1):
            cnt3+= 1
            
        N_corr_arr = np.array(N_corr_arr)
        mean_N,med_N, std_N = sigma_clipped_stats(N_corr_arr, cenfunc=np.median)
        
             
        #Now correct for PSF and use monte carlo if needed
        corr_sigxx_arr=[]
        corr_sigyy_arr=[]
        corr_sigxy_arr=[]
        wtArr = np.array(wtArr)
        wtArr = np.ones(len(wtArr))/np.sum(wtArr)
        #print (wtArr)
        wtArr = wtArr**2#(Inverse variance weight)
        wtArr = wtArr/np.sum(wtArr)
        #print (list(wtArr))
        
        
        median_size= np.nanmedian(size_sf_arr)
        sigxx_err_arr = np.array(sigxx_err_arr)
        sigyy_err_arr = np.array(sigyy_err_arr)
        sigxy_err_arr = np.array(sigxy_err_arr)
        sigxx_arr = np.array(sigxx_arr)
        sigyy_arr = np.array(sigyy_arr)
        sigxy_arr = np.array(sigxy_arr)
        sigxx_err_arr = np.array(sigxx_err_arr)
        sigyy_err_arr = np.array(sigyy_err_arr)
        sigxy_err_arr = np.array(sigxy_err_arr)
        psfsigxx_arr = np.array(psfsigxx_arr)
        psfsigyy_arr = np.array(psfsigyy_arr)
        psfsigxy_arr = np.array(psfsigxy_arr)
        
        
        mean_xx, med_xx, std_xx = sigma_clipped_stats(sigxx_arr - psfsigxx_arr, cenfunc=np.median)
        mean_yy, med_yy, std_yy = sigma_clipped_stats(sigyy_arr - psfsigyy_arr, cenfunc=np.median)
        mean_xy, med_xy, std_xy = sigma_clipped_stats(sigxy_arr - psfsigxy_arr, cenfunc=np.median)
        
        validIndices = np.where((wtArr >0) & (N_corr_arr< (med_N+3*std_N)) & (N_corr_arr> (med_N-3*std_N)) 
                                & (N_corr_arr> 0 ) )[0]
        
        for k in range(len(sigxx_arr)):
            corr_xx = sigxx_arr[k] - psfsigxx_arr[k]
            corr_yy = sigyy_arr[k] - psfsigyy_arr[k]
            corr_xy = sigxy_arr[k] - psfsigxy_arr[k]
            corr_sigxx_arr.append(corr_xx)
            corr_sigyy_arr.append(corr_yy)
            corr_sigxy_arr.append(corr_xy)
            if(k not in validIndices and len(wtArr)>10):
                wtArr[k] = 0
           
        print (wtArr)
        wtArr = wtArr/np.sum(wtArr)
        corr_sigxx_arr = np.array(corr_sigxx_arr)
        corr_sigyy_arr = np.array(corr_sigyy_arr)
        corr_sigxy_arr = np.array(corr_sigxy_arr)
        

                
            
        
        
        
        xx_coadd = ir_coadd_data [j,35]-ir_coadd_data [j,38]
        yy_coadd = ir_coadd_data [j,36]-ir_coadd_data [j,39]
        xy_coadd = ir_coadd_data [j,37]+ir_coadd_data [j,40]
        e1_coadd = (xx_coadd-yy_coadd)/(xx_coadd+yy_coadd)
        e2_coadd = 2*xy_coadd/(xx_coadd+yy_coadd)
        
        #print (e1, e1_coadd, j)
        
        #print (corr_sigxx-ir_coadd_data [j,35]+ir_coadd_data [j,38], corr_sigyy-ir_coadd_data [j,36]+ir_coadd_data [j,39])
        
        master_frame[j, 0, 0:len(corr_sigxx_arr)] = np.ones(len(corr_sigxx_arr)) * ir_coadd_data[j,10]
        master_frame[j, 1, 0:len(corr_sigxx_arr)] = np.ones(len(corr_sigxx_arr)) *ir_coadd_data[j,11]
        master_frame[j, 2, 0:len(corr_sigxx_arr)] = (corr_sigxx_arr - corr_sigyy_arr)/(corr_sigxx_arr + corr_sigyy_arr)
        master_frame[j, 3, 0:len(corr_sigxx_arr)] = 2*corr_sigxy_arr/(corr_sigxx_arr + corr_sigyy_arr)
        master_frame[j, 4, 0:len(corr_sigxx_arr)] = np.ones(len(corr_sigxx_arr)) *redShiftArr[j]
        master_frame[j,5, 0] = len(validIndices)
        master_frame[j,6,0:len(corr_sigxx_arr)] = corr_sigxx_arr
        master_frame[j,7,0:len(corr_sigxx_arr)] = corr_sigyy_arr
        master_frame[j,8,0:len(corr_sigxx_arr)] = corr_sigxy_arr
        master_frame[j,9, 0:len(corr_sigxx_arr)] = np.ones(len(corr_sigxx_arr)) *ir_coadd_data [j,3]
        master_frame[j,10,0:len(corr_sigxx_arr)] = sigxx_err_arr
        master_frame[j,11,0:len(corr_sigxx_arr)] = sigyy_err_arr
        master_frame[j,12,0:len(corr_sigxx_arr)] = sigxy_err_arr
        master_frame[j,13,0:len(corr_sigxx_arr)] = np.sqrt( sigxx_err_arr**2 * 4* (corr_sigxx_arr**2+corr_sigyy_arr**2)/(corr_sigxx_arr+corr_sigyy_arr)**4)
        master_frame[j,14,0:len(corr_sigxx_arr)] = np.sqrt( (sigxx_err_arr**2 * 8* (corr_sigxy_arr**2)/(corr_sigxx_arr+corr_sigyy_arr)**4) + (sigxy_err_arr**2 * 4/(corr_sigxx_arr+corr_sigyy_arr)**2))
        #master_frame[j,15,0:len(corr_sigxx_arr)] = np.sqrt((e1_coadd**2 *master_frame[j,13, 0:len(corr_sigxx_arr)]**2 / (e1_coadd**2 + e2_coadd**2)) + (e2_coadd**2 *master_frame[j,14, 0:len(corr_sigxx_arr)]**2 / (e1_coadd**2 + e2_coadd**2)) )
        
        
        master_frame[j,15,0:len(corr_sigxx_arr)] = np.sqrt((master_frame[j, 2, 0:len(corr_sigxx_arr)]**2 *master_frame[j,13, 0:len(corr_sigxx_arr)]**2 / (master_frame[j, 2, 0:len(corr_sigxx_arr)]**2 + master_frame[j, 3, 0:len(corr_sigxx_arr)]**2)) + (master_frame[j, 3, 0:len(corr_sigxx_arr)]**2 *master_frame[j,14, 0:len(corr_sigxx_arr)]**2 / (master_frame[j, 2, 0:len(corr_sigxx_arr)]**2 + master_frame[j, 3, 0:len(corr_sigxx_arr)]**2)) )

        

     
            
    np.save(outFile, master_frame)
    
    
    
    
    
    
from astropy.io import fits
import numpy as np
import subprocess, os, shutil
band ='ir'
img = '/scratch/bell/dutta26/wiyn_sim/wted_coadds/ir_coadd_wt.fits'
med_img = '/scratch/bell/dutta26/wiyn_sim/median_coadds/ir_coadd_med_crop.fits'
sextractorFile = '/scratch/bell/dutta26/wiyn_sim/coadd_sim.cat'
starGal_out_file = '/home/dutta26/codes/wiyn_wl_sim/source_list.pk1'
coadd_detect_file = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_' + str(band)+'.npy'
ebModeFile = '/scratch/bell/dutta26/wiyn_sim/master_arr_coaddMC_sfMC_all.npy'
plotLoc = '/scratch/bell/dutta26/wiyn_sim/plot/phosim/'
bandLoc = '/scratch/bell/dutta26/wiyn_sim/'
ir_band_data = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy'



EandB_all( '/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy',  '/home/dutta26/codes/wiyn_wl_sim/coaddSc_r.npy', '/home/dutta26/codes/wiyn_wl_sim/coaddSc_i.npy', None, ebModeFile, img, '/scratch/bell/dutta26/wiyn_sim/r_withMC.npy', '/scratch/bell/dutta26/wiyn_sim/i_withMC.npy')

