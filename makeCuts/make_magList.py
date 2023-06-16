#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 09:14:42 2022

@author: dutta26
"""

import numpy as np
import pandas as pd
import helper,sys,os
from astropy.io import fits

source_df = pd.read_pickle('/home/dutta26/codes/source_list.pk1')
raList = np.array(source_df['ra'])
decList = np.array(source_df['dec'])


coadd_data_g = pd.read_pickle('/home/dutta26/codes/coaddSc_g_.pk1')
coadd_data_g = np.array(coadd_data_g)

coadd_data_r = pd.read_pickle('/home/dutta26/codes/coaddSc_r_.pk1')
coadd_data_r = np.array(coadd_data_r)

coadd_data_i = pd.read_pickle('/home/dutta26/codes/coaddSc_i_.pk1')
coadd_data_i = np.array(coadd_data_i)

coadd_data_z = pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1')
coadd_data_z = np.array(coadd_data_z)

coadd_data_u = pd.read_pickle('/home/dutta26/codes/coaddSc_u_.pk1')
coadd_data_u = np.array(coadd_data_u)

a,b = np.shape(coadd_data_g)

band_coadd_df = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
band_coadd_xx = np.array(band_coadd_df['xx'])
band_coadd_yy = np.array(band_coadd_df['yy'])




tempInd = np.where((coadd_data_u[:, 3] <= 0))
xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/abell_u_coadd_withBkg.fits')
f=fits.open('/scratch/bell/dutta26/abell_2390/abell_u_coadd_withBkg.fits')
data = np.array(f[0].data)
f.close()

mag_u = np.ones(len(coadd_data_u)) *-99
err_u = np.ones(len(coadd_data_u)) *-99
for j in range(len(coadd_data_u)):
    if(j in tempInd[0]):
        continue
    x = int(xList[j])
    y = int(yList[j])
    B = np.mean(data[y-50:y+50, x-50:x+50])
    A = 2*np.pi*np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
    error_in_flux = np.sqrt(4*A*B +coadd_data_u[j, 3] )
    if(error_in_flux == None or np.isnan(error_in_flux)):
        print (j)
        break
    flux =   coadd_data_u[j, 3]  *np.sqrt(len(os.listdir('/scratch/bell/dutta26/abell_2390/u/')) / 2)                    
    mag_u[j] = 25 - 2.5*np.log10(coadd_data_u[j, 3]/60) - 4.445378
    err_u[j] = (2.5/np.log(10)) * (error_in_flux/flux  )
    
mag_u[tempInd[0]] = -99
err_u[tempInd[0]] = -99
print (len(tempInd[0]))







tempInd = np.where((coadd_data_g[:, 3] <= 0))
xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/abell_g_coadd_withBkg.fits')
f=fits.open('/scratch/bell/dutta26/abell_2390/abell_g_coadd_withBkg.fits')
data = np.array(f[0].data)
f.close()

mag_g = np.ones(len(coadd_data_g)) *-99
err_g = np.ones(len(coadd_data_g)) *-99
for j in range(len(coadd_data_g)):
    if(j in tempInd[0]):
        continue
    
    x = int(xList[j])
    y = int(yList[j])
    B = np.mean(data[y-50:y+50, x-50:x+50])
    A = 2*np.pi*np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
    error_in_flux = np.sqrt(4*A*B +coadd_data_g[j, 3] )
    if(error_in_flux == None or np.isnan(error_in_flux)):
        print (j)
        break
    flux =   coadd_data_g[j, 3]  *np.sqrt(len(os.listdir('/scratch/bell/dutta26/abell_2390/g/')) / 2) 
    
    mag_g[j] = 25 - 2.5*np.log10(coadd_data_g[j, 3]/60) - 4.445378
    err_g[j] = (2.5/np.log(10)) * (error_in_flux/flux)
    
mag_g[tempInd[0]] = -99
err_g[tempInd[0]] = -99
print (len(tempInd[0]))






tempInd = np.where((coadd_data_r[:, 3] <= 0))
xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/abell_r_coadd_withBkg.fits')
f=fits.open('/scratch/bell/dutta26/abell_2390/abell_r_coadd_withBkg.fits')
data = np.array(f[0].data)
f.close()

mag_r = np.ones(len(coadd_data_r)) *-99
err_r = np.ones(len(coadd_data_r)) *-99
for j in range(len(coadd_data_r)):
    if(j in tempInd[0]):
        continue
    x = int(xList[j])
    y = int(yList[j])
    B = np.mean(data[y-50:y+50, x-50:x+50])
    A = 2*np.pi*np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
    error_in_flux = np.sqrt(4*A*B +coadd_data_r[j, 3] )
    if(error_in_flux == None or np.isnan(error_in_flux)):
        print (j)
        break
    flux =   coadd_data_r[j, 3]  *np.sqrt(len(os.listdir('/scratch/bell/dutta26/abell_2390/r/')) / 2) 
    mag_r[j] = 25 - 2.5*np.log10(coadd_data_r[j, 3]/60) - 4.445378
    err_r[j] = (2.5/np.log(10)) * (error_in_flux/flux)
    
mag_r[tempInd[0]] = -99
err_r[tempInd[0]] = -99
print (len(tempInd[0]))







tempInd = np.where((coadd_data_i[:, 3] <= 0))
xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/abell_i_coadd_withBkg.fits')
f=fits.open('/scratch/bell/dutta26/abell_2390/abell_i_coadd_withBkg.fits')
data = np.array(f[0].data)
f.close()
mag_i = np.ones(len(coadd_data_i)) *-99
err_i = np.ones(len(coadd_data_i)) *-99
for j in range(len(coadd_data_i)):
    if(j in tempInd[0]):
        continue
    x = int(xList[j])
    y = int(yList[j])
    B = np.mean(data[y-50:y+50, x-50:x+50])
    A = 2*np.pi*np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
    error_in_flux = np.sqrt(4*A*B +coadd_data_i[j, 3] )
    if(error_in_flux == None or np.isnan(error_in_flux)):
        print (j)
        break
    flux =   coadd_data_i[j, 3]  *np.sqrt(len(os.listdir('/scratch/bell/dutta26/abell_2390/i/')) / 2)   
    mag_i[j] = 25 - 2.5*np.log10(coadd_data_i[j, 3]/60) - 4.445378
    err_i[j] = (2.5/np.log(10)) * (error_in_flux/flux)
    
mag_i[tempInd[0]] = -99
err_i[tempInd[0]] = -99
print (len(tempInd[0]))





tempInd = np.where((coadd_data_z[:, 3] <= 0))
xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/abell_z_coadd_withBkg.fits')
f=fits.open('/scratch/bell/dutta26/abell_2390/abell_z_coadd_withBkg.fits')
data = np.array(f[0].data)
f.close()
mag_z = np.ones(len(coadd_data_z)) *-99
err_z = np.ones(len(coadd_data_z)) *-99
for j in range(len(coadd_data_z)):
    if(j in tempInd[0]):
        continue
    x = int(xList[j])
    y = int(yList[j])
    B = np.mean(data[y-50:y+50, x-50:x+50])
    A = 2*np.pi*np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
    error_in_flux = np.sqrt(4*A*B +coadd_data_z[j, 3] )
    if(error_in_flux == None or np.isnan(error_in_flux)):
        print (j)
        break
    flux =   coadd_data_z[j, 3]  * np.sqrt(len(os.listdir('/scratch/bell/dutta26/abell_2390/z/')) / 2) 
    mag_z[j] = 25 - 2.5*np.log10(coadd_data_z[j, 3]/60) - 4.445378
    err_z[j] = (2.5/np.log(10)) * (error_in_flux/flux)
    
mag_z[tempInd[0]] = -99
err_z[tempInd[0]] = -99
print (len(tempInd[0]))



















# =============================================================================
# tempInd = np.where((coadd_data_g[:, 18] == 0) | (coadd_data_g[:, 18] >32 ) | (coadd_data_g[:, 18] <10 ))
# mag_g = coadd_data_g[:, 18]
# mag_g[tempInd[0]] = -99
# err_g = np.sqrt(np.power(10,0.4*(20.554622-coadd_data_g[:, 18])))/(10*np.power(10,0.4*(20.554622-coadd_data_g[:, 18])))
# err_g[tempInd[0]] = -99
# print (len(tempInd[0]))
# 
# tempInd = np.where((coadd_data_r[:, 18] == 0) | (coadd_data_r[:, 18] >32 ) | (coadd_data_r[:, 18] <10 ))
# mag_r = coadd_data_r[:, 18]
# mag_r[tempInd[0]] = -99
# err_r = np.sqrt(np.power(10,0.4*(20.554622-coadd_data_r[:, 18])))/(10*np.power(10,0.4*(20.554622-coadd_data_r[:, 18])))
# err_r[tempInd[0]] = -99
# print (len(tempInd[0]))
# 
# 
# tempInd = np.where((coadd_data_i[:, 18] == 0) | (coadd_data_i[:, 18] >32 ) | (coadd_data_i[:, 18] <10 ))
# mag_i = coadd_data_i[:, 18]
# mag_i[tempInd[0]] = -99
# err_i = np.sqrt(np.power(10,0.4*(20.554622-coadd_data_i[:, 18])))/(10*np.power(10,0.4*(20.554622-coadd_data_i[:, 18])))
# err_i[tempInd[0]] = -99
# print (len(tempInd[0]))
# 
# 
# tempInd = np.where((coadd_data_z[:, 18] == 0) | (coadd_data_z[:, 18] >32 ) | (coadd_data_z[:, 18] <10 ))
# mag_z = coadd_data_z[:, 18]
# mag_z[tempInd[0]] = -99
# err_z = np.sqrt(np.power(10,0.4*(20.554622-coadd_data_z[:, 18])))/(10*np.power(10,0.4*(20.554622-coadd_data_z[:, 18])))
# err_z[tempInd[0]] = -99
# print (len(tempInd[0]))
# =============================================================================


outfile = '/home/dutta26/A2390_mag1.in'
f = open(outfile, mode='w+')

for j in range(a):
    f.write(str(j) + ' ' + str(np.round(mag_u[j], 3))+ ' '+ str(np.round(err_u[j], 3))+ ' '+
            str(np.round(mag_g[j], 3))+ ' '+ str(np.round(err_g[j],3))+ ' ' + 
            str(np.round(mag_r[j],3))+ ' '+ str(np.round(err_r[j], 3))+ ' '+
            str(np.round(mag_i[j],3))+ ' '+ str(np.round(err_i[j],3))+ ' '+
            str(np.round(mag_z[j],3))+ ' '+ str(np.round(err_z[j],3))+ ' '+' \n')
    
f.close()


