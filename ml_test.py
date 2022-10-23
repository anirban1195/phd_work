#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 08:38:09 2021

@author: dutta26
"""
import numpy as np
import pandas as pd
band = 'r'
frame_data = np.load('/home/dutta26/codes/singleFrame_'+band+'.npy')
coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')

ir_coadd_data = np.array(ir_coadd_data)
coadd_data = np.array(coadd_data)

#Find all the stars 
indexList = np.where(coadd_data[:,2] == 1)[0]
percent_flux =[]

airmass =[]
pos =[]
size =[]
see =[]
delPos = []
mjd = []
bkg = []
zp = []
wt=[]
ind =[]
for j in indexList[11:200]:
    coadd_flux = coadd_data[j, 3]
    
    a = frame_data[:,j, :]
    temp = np.where( (a[:,15] ==0) & (a[:,16] ==0) & (a[:,17] ==0) &(a[:,2] ==1) )[0]
    a = a[ np.where( (a[:,15] ==0) & (a[:,16] ==0) & (a[:,17] ==0) &(a[:,2] ==1) ), :]
    a=a[0,:,:]
    
    p,q = np.shape(a)
    for k in range(p):
        percent_flux.append( ((a[k,3] - coadd_flux)*100) / coadd_flux)
        ind.append(temp[k])
        airmass.append( a[k,23]*10 )
        pos.append( np.sqrt(a[k,13]**2 + a[k,14]**2)/100 )
        size.append( np.sqrt( a[k,10]**2 + a[k,11]**2 ) )
        see.append(a[k,19] * 10)
        delPos.append( np.sqrt(a[k,4]**2 + a[k,5]**2 ) *10 )
        mjd.append((a[k,22] - 50000)/100)
        bkg.append(a[k,18]/ 10)
        zp.append(a[k, 20] -20)
        wt.append( np.power(10,(a[k, 20]-25)/2.5)/((np.sqrt(a[k,18])*a[k,19])**2))
    break
    