#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 12:29:11 2023

@author: dutta26
"""

import numpy as np
import pandas as pd
import helper,sys,os
from astropy.io import fits


def getMagErr(band, sourceFile ):
    source_df = pd.read_pickle(sourceFile)
    raList = np.array(source_df['ra'])
    decList = np.array(source_df['dec'])
    coadd_data_band =  np.load('/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.npy') #####
    coadd_data_ir = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')######
    a,b = np.shape(coadd_data_ir)
    
    fileList=[]
    f=open('/home/dutta26/codes/makeWeights/fileList_swarp_'+band+'.ascii')#####
    content = f.readlines()
    f.close()
    for j in range(len(content)):
        fileList.append(content[j][:-1])
    
    totBkg = 0
    totScale= 0
    for files in fileList:
        if('temp' in files or 'weight' in files):
            continue
        f=fits.open(files)
        back_h = float((f[0].header)['SKYBG'])
        flux_scale = float((f[0].header)['FLXSCALE']) 
        f.close()
        totBkg += back_h
        print (flux_scale)
        totScale += 1/flux_scale
    
    
    mag = np.ones(len(coadd_data_band)) *-99
    err = np.ones(len(coadd_data_band)) *-99
    for j in range(len(coadd_data_band)):
        #print ('aa')
        flux =   coadd_data_band[j, 3] *totScale
        size = np.sqrt(coadd_data_ir[j,35] + coadd_data_ir[j,36])
        error_in_flux = np.sqrt(flux + 4*np.pi*size**2 * totBkg )
        
        if(error_in_flux == None or np.isnan(error_in_flux) or flux <= 0 or flux==None or np.isnan(flux) or np.isinf(flux) or np.isinf(error_in_flux)):
            #print (j)
            continue
                      
        mag[j] = 25 - 2.5*np.log10(coadd_data_band[j, 3]) 
        err[j] = (2.5/np.log(10)) * (error_in_flux/flux  )
        if(err[j]< 0.01):
            err[j]=0.01
        
    return mag,err





sourceFile = '/home/dutta26/codes/source_list.pk1'
coadd_data_ir = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
a,b = np.shape(coadd_data_ir)
bandList = ['u', 'g', 'r', 'i', 'z']
mag_arr = np.ones((10, a), dtype = np.float32) * -99
cnt = 0
for band in bandList:
    mag, err= getMagErr(band, sourceFile)
    mag_arr[cnt,:], mag_arr[cnt+ 1, :]= mag, err
    cnt += 2
    


outfile = '/home/dutta26/A2390_mag4.in'
f = open(outfile, mode='w+')

for j in range(a):
    f.write(str(j) + ' ' + str(np.round(mag_arr[0,j], 3))+ ' '+ str(np.round(mag_arr[1,j], 3))+ ' '+
            str(np.round(mag_arr[2,j], 3))+ ' '+ str(np.round(mag_arr[3,j],3))+ ' ' + 
            str(np.round(mag_arr[4,j],3))+ ' '+ str(np.round(mag_arr[5,j], 3))+ ' '+
            str(np.round(mag_arr[6,j],3))+ ' '+ str(np.round(mag_arr[7,j],3))+ ' '+
            str(np.round(mag_arr[8,j],3))+ ' '+ str(np.round(mag_arr[9,j],3))+ ' '+' \n')
    
f.close()


