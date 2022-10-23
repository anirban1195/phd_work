#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 12:38:41 2022

@author: dutta26
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess
from scipy.special import erf

sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test_10star_1.npy')
coadd_df = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')

loc = np.where((coadd_df[:,3]>1000) & (coadd_df[:,3]<1100) & (coadd_df[:,2] == 0 ) )[0]
e1Arr=[]
e2Arr=[]
e1_devArr=[]
e2_devArr =[]
cnt1 = cnt2= cnt3= cnt4 = cnt5= cnt6 = cnt7 =cnt8 = cnt9 = 0
a1=a2=a3=a4=a5=a6=0
indexArr=[]
for j in loc:
    corr_xx = coadd_df[j,7] - coadd_df[j,38]
    corr_yy = coadd_df[j,8] - coadd_df[j,39]
    corr_xy = coadd_df[j,9] - coadd_df[j,40]
    if(corr_xx< 0 or corr_yy<0):       
        continue 
   
           
    e1 = (corr_xx - corr_yy)/(corr_xx+corr_yy)
    e2 = 2*corr_xy/(corr_xx+corr_yy)
    
    if(np.abs(e1)>1 or np.abs(e2)>1):
        a6+= 1
        continue
    #sys.exit()
    e1Arr.append(e1)
    e2Arr.append(e2)
    
    temp1 =[]
    temp2 =[]
    tempwt =[]
    indexArr=[]
    for k in range(80):
        #Makre sure they have been measured
        if(sf_df[k,j,12]!= 99 and sf_df[k,j,12]!= 1 ):
            continue
        if(sf_df[k,j,38] == -99 ):
            continue
        
        
        if(sf_df[k,j,12]== 99):
            if(sf_df[k,j,7]<=0 or sf_df[k,j,8]<=0):
                a1+= 1
                continue
            if(np.abs(sf_df[k,j,7])>=30 or np.abs(sf_df[k,j,8])>=30):
                a2+= 1
                continue
            
        if(sf_df[k,j,12]== 1):
            if(sf_df[k,j,35]<=0 or sf_df[k,j,36]<=0):
                a3 += 1
                continue
            if(np.abs(sf_df[k,j,35])>= 30 or np.abs(sf_df[k,j,36])>=30):
                a4 += 1
                continue
            
        
        #print (k,sf_df[k,j,12])
        if(sf_df[k,j,12]== 99):
            m_xx = sf_df[k,j,7]
            m_yy = sf_df[k,j,8]
            m_xy = sf_df[k,j,9]
            B = sf_df[k,j,6]
            N = sf_df[k,j,3]
            cnt2 += 1
            
        elif(sf_df[k,j,12]== 1):
            m_xx = sf_df[k,j,35]
            m_yy = sf_df[k,j,36]
            m_xy = sf_df[k,j,37]
            B = sf_df[k,j,34]
            N = sf_df[k,j,31]
            cnt3 += 1
        if(N<= 0):
            continue 
        
        
        cnt1 += 1
        corr_xx = m_xx - sf_df[k,j,38]
        corr_yy = m_yy - sf_df[k,j,39]
        corr_xy = m_xy - sf_df[k,j,40]
        temp = corr_xx +corr_yy - 2*np.abs(corr_xy)
        snr = (N/ np.sqrt( N + 4*B*np.sqrt(m_xx*m_yy)*3.14 ))
        if(np.isnan(snr)):
            sys.exit()
            
        if(corr_xx< 0 or corr_yy<0 or temp<0):
            if(sf_df[k,j,12]== 99):
                cnt4 += 1
            elif(sf_df[k,j,12]== 1):
                cnt5+= 1
            temp_xx = coadd_df[j,7] - coadd_df[j,38] + sf_df[k,j,38]
            temp_yy = coadd_df[j,8] - coadd_df[j,39] + sf_df[k,j,39]
            temp_xy = coadd_df[j,9] - coadd_df[j,40] + sf_df[k,j,40]
            A = np.pi * np.sqrt(2*temp_xx* 2*temp_yy - (2*temp_xy)**2)
            s = np.sqrt( (A/(np.pi*N) + 4*A**2 * B/(np.pi * N**2)) ) * np.sqrt(temp_xx* temp_yy - (temp_xy)**2)
            #print (s, corr_xx, corr_yy,coadd_df[j,41], '*****' )
            err1 = 1.16*(2*temp_xx)/np.sqrt(N)
            err2 = 1.5 *(2*temp_xx) * np.sqrt(4*A*B)/ N
            s_xx = (1+e1)*s
            
            err1 = 1.16*(2*temp_yy)/np.sqrt(N)
            err2 = 1.5 *(2*temp_yy) * np.sqrt(4*A*B)/ N
            s_yy = (1-e1)*s
            
            err1 = 0.8*(2*np.sqrt(temp_yy) * np.sqrt(temp_xx))/np.sqrt(N)
            err2 = 1.5 *(2*np.sqrt(temp_yy) * np.sqrt(temp_xx)) * np.sqrt(4*A*B)/ N
            s_xy = abs(e2)*s
            
            corr_xx, corr_yy, corr_xy , success = helper.correct(m_xx, m_yy, m_xy, sf_df[k,j,38], sf_df[k,j,39], sf_df[k,j,40], s_xx,s_yy,s_xy, 
                                                           sf_df[k,j,41], sf_df[k,j,42], sf_df[k,j,43])
            #print (success, corr_xx, corr_yy)
            
            if(success<25 ):
                print (s, corr_xx, m_xx - sf_df[k,j,38],m_xx ,sf_df[k,j,38], j,k)
                pass
            else:
                if(sf_df[k,j,12]== 99):
                    cnt6 += 1
                elif(sf_df[k,j,12]== 1):
                    cnt7+= 1
                e1_sf1 = (corr_xx - corr_yy)/(corr_xx+corr_yy)
                e2_sf1 = 2*corr_xy/(corr_xx+corr_yy)
                temp1.append(e1_sf1)
                temp2.append(e2_sf1)
                indexArr.append(k)
                if(sf_df[k,j,12]== 99):
                    tempwt.append(snr)
                else:
                    tempwt.append(snr/5)
                
        else:
            if(sf_df[k,j,12]== 99):
                
                cnt8 += 1
            elif(sf_df[k,j,12]== 1):
                cnt9+= 1
            e1_sf1 = (corr_xx - corr_yy)/(corr_xx+corr_yy)
            e2_sf1 = 2*corr_xy/(corr_xx+corr_yy)
            temp1.append(e1_sf1)
            temp2.append(e2_sf1)
            indexArr.append(k)
            if(sf_df[k,j,12]== 99):
                tempwt.append(snr)
            else:
                tempwt.append(snr/5)
        
    temp1 = np.array(temp1)    
    temp2 = np.array(temp2)
    tempwt = np.array(tempwt)
    mean_e1= np.sum(temp1*tempwt)/np.sum(tempwt)
    mean_e2= np.sum(temp2*tempwt)/np.sum(tempwt)
    e1_devArr.append(e1-mean_e1)
    e2_devArr.append(e2-mean_e2)
    if(abs(e1-mean_e1)> 0.2):
        break


    
n, bins, patches = plt.hist(x=e1_devArr, bins='auto', color='r',
                               alpha=0.5, rwidth=0.85, density=True)    
    
    
n, bins, patches = plt.hist(x=e2_devArr, bins='auto', color='b',
                               alpha=0.5, rwidth=0.85, density=True) 
       
    
plt.title("SNR ~ 7 ")    
print (sigma_clipped_stats(e1_devArr), sigma_clipped_stats(e2_devArr))    
print (len(loc), len(e1_devArr)/80, len(e2_devArr)/80)