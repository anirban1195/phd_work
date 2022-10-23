#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 13:34:08 2022

@author: dutta26
"""


from astroquery.sdss import SDSS

from astropy import coordinates as coords
from astroquery.ned import Ned
from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy import coordinates
import numpy as np
import pandas as pd

SDSS.TIMEOUT=500

coadd_data_g = pd.read_pickle('/home/dutta26/codes/coaddSc_new_g_.pk1')
coadd_data_g = np.array(coadd_data_g)

coadd_data_r = pd.read_pickle('/home/dutta26/codes/coaddSc_new_r_.pk1')
coadd_data_r = np.array(coadd_data_r)

coadd_data_i = pd.read_pickle('/home/dutta26/codes/coaddSc_new_i_.pk1')
coadd_data_i = np.array(coadd_data_i)

coadd_data_z = pd.read_pickle('/home/dutta26/codes/coaddSc_new_z_.pk1')
coadd_data_z = np.array(coadd_data_z)

coadd_data_u = pd.read_pickle('/home/dutta26/codes/coaddSc_new_u_.pk1')
coadd_data_u = np.array(coadd_data_u)

a,b = np.shape(coadd_data_g)



tempInd = np.where((coadd_data_u[:, 18] == 0) | (coadd_data_u[:, 18] >32 ) | (coadd_data_u[:, 18] <10 ))
mag_u = coadd_data_u[:, 18]
mag_u[tempInd[0]] = -99

print (len(tempInd[0]))

tempInd = np.where((coadd_data_g[:, 18] == 0) | (coadd_data_g[:, 18] >32 ) | (coadd_data_g[:, 18] <10 ))
mag_g = coadd_data_g[:, 18]
mag_g[tempInd[0]] = -99

print (len(tempInd[0]))

tempInd = np.where((coadd_data_r[:, 18] == 0) | (coadd_data_r[:, 18] >32 ) | (coadd_data_r[:, 18] <10 ))
mag_r = coadd_data_r[:, 18]
mag_r[tempInd[0]] = -99

print (len(tempInd[0]))

tempInd = np.where((coadd_data_z[:, 18] == 0) | (coadd_data_z[:, 18] >32 ) | (coadd_data_z[:, 18] <10 ))
mag_z = coadd_data_z[:, 18]
mag_z[tempInd[0]] = -99

print (len(tempInd[0]))


tempInd = np.where((coadd_data_i[:, 18] == 0) | (coadd_data_i[:, 18] >32 ) | (coadd_data_i[:, 18] <10 ))
mag_i = coadd_data_i[:, 18]
mag_i[tempInd[0]] = -99

print (len(tempInd[0]))


uDiff=[]
gDiff=[]
rDiff=[]
iDiff=[]
zDiff=[]
for j in range(a):
    ra = coadd_data_r[j,0]
    dec = coadd_data_r[j,1]
    if(ra == 0 or dec ==0):
        continue
    print (j)
    if( mag_g[j] >24 or mag_r[j]>24 or mag_i[j] >24 ):
        continue
    
    
    pos = coords.SkyCoord(ra=ra, dec=dec,
                          unit=(u.deg, u.deg), frame='icrs')
    result_table = SDSS.query_region(pos, radius=0.0001 * u.deg, photoobj_fields=['ra', 'dec', 'cmodelMag_u', 'cmodelMag_g', 'cmodelMag_r', 'cmodelMag_i', 'cmodelMag_z'])
    if(result_table == None):
        continue
    
    

    df = result_table.to_pandas()
    df = df.sort_values(by=['cmodelMag_g'])
    df = df.iloc[0]
    
    
    if(mag_u[j] != -99 and df['cmodelMag_u'] != -9999):
        uDiff.append(mag_u[j] - df['cmodelMag_u'])
        
    
    if(mag_g[j] != -99 and df['cmodelMag_g'] != -9999):
        gDiff.append(mag_g[j] - df['cmodelMag_g'])  
        
    if(mag_r[j] != -99 and df['cmodelMag_r'] != -9999):
        rDiff.append(mag_r[j] - df['cmodelMag_r'])
        
    if(mag_i[j] != -99 and df['cmodelMag_i'] != -9999):
        iDiff.append(mag_i[j] - df['cmodelMag_i'])
        
    if(mag_z[j] != -99 and df['cmodelMag_z'] != -9999):
        zDiff.append(mag_z[j] - df['cmodelMag_z'])
        
uDiff = np.array(uDiff)
gDiff = np.array(gDiff)
rDiff = np.array(rDiff)
iDiff = np.array(iDiff)
zDiff = np.array(zDiff)          
np.save('/scratch/halstead/d/dutta26/abell_2390/diff_u.npy', uDiff)
np.save('/scratch/halstead/d/dutta26/abell_2390/diff_g.npy', gDiff)
np.save('/scratch/halstead/d/dutta26/abell_2390/diff_r.npy', rDiff)
np.save('/scratch/halstead/d/dutta26/abell_2390/diff_i.npy', iDiff)
np.save('/scratch/halstead/d/dutta26/abell_2390/diff_z.npy', zDiff)
