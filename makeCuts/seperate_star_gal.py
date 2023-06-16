#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 22:11:19 2021

@author: dutta26
"""


import astropy.units as u
import helper
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from astroquery.gaia import Gaia
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import sys

Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
# =============================================================================
# catalog = str(sys.argv[1])
# outFile = str(sys.argv[2])
# imgFile = str(sys.argv[3])
# =============================================================================
catalog = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/abell_2390_sextract.cat'
outFile = '/home/dutta26/codes/source_list.pk1'
imgFile = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd_crop.fits'

f=fits.open(imgFile)
img = np.array(f[0].data)
sizey, sizex = np.shape(img)

#Find the max, min and mid Ra Dec
w=wcs.WCS(f[0].header)
f.close()
[[startRa, startDec]]  = w.wcs_pix2world([[0,0]], 0)
[[endRa, endDec]] = w.wcs_pix2world([[sizex, sizey]], 0)
midRa = 0.5*(startRa +endRa )
midDec = 0.5* (startDec + endDec)

#Get stars from GAIA within a degree of mid ra dec
coord = SkyCoord(ra=midRa , dec=midDec, unit=(u.degree, u.degree), frame='icrs')
width = u.Quantity(np.abs(startRa - endRa ), u.deg)
height = u.Quantity(np.abs(startDec - endDec), u.deg)
Gaia.ROW_LIMIT = -1
r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
r.pprint()

df = r.to_pandas()

#sys.exit()

#df.drop(df[df['parallax_over_error'] < 0.5].index, inplace = True)
df.drop(df[df['phot_g_mean_mag'] < 14].index, inplace = True)

#df1 = df[['ra', 'dec', 'parallax_over_error']]
df2 = df[df['parallax_over_error'].notna()]

#Read the catalog 
with open(catalog) as f:
    content = f.readlines()

Raindex=5
Decindex=6
xxIndex = 7
yyIndex = 8
xyIndex = 9

#Make an array containing x and y indices 
raList =[]
decList =[]
sex_xx_list = []
sex_yy_list = []
sex_xy_list = []
star_list =[]
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    
    raList.append(float(content[j].split()[Raindex])) 
    decList.append(float(content[j].split()[Decindex])) 
    sex_xx_list.append(float(content[j].split()[xxIndex]))
    sex_yy_list.append(float(content[j].split()[yyIndex]))
    sex_xy_list.append(float(content[j].split()[xyIndex]))
    star_list.append(0)
    
    
#Make a pandas dataframe for the sources . 1 for star and 0 for galaxy
a=np.array([raList, decList, sex_xx_list, sex_yy_list, sex_xy_list, star_list])
b=np.swapaxes(a, 0,1)
df_source = pd.DataFrame(b,  columns = ['ra', 'dec', 'sex_xx', 'sex_yy', 'sex_xy', 'star_bool'])

    




finalMatchIndexArr =[]
for j in range(len(df2)):
    
#for j in [2947, 2954]:
    ra = float( df2['ra'][j:j+1])
    dec = float( df2['dec'][j:j+1])
    
    fluxMatch=0
    matchFlag = 0 
    matchIndexArr=[]
    distArr =[]
    #print ('$$$$$$$$$$$$$$$$$$')
    for k in range(len(raList)):
        if(np.abs(raList[k] - ra) < 0.003 and np.abs(decList[k] - dec) < 0.003):
            matchIndexArr.append(k)
            distArr.append(np.sqrt((raList[k] - ra)**2 + (decList[k] - dec)**2))
    
    
    #If no match is found       
    if(len(matchIndexArr) == 0):
        continue
        
    #If only 1 match is found and is close 
    if(len(matchIndexArr) == 1 and distArr[0]<= 0.0002):
        #print ('aaaaaaa')
        df_source['star_bool'][matchIndexArr[0]] = 1
        continue
    #If multiple matches are found 
    temp_dist = np.sort(distArr)
    #print (temp_dist)
    #if(temp_dist[0]< 0.0002 and temp_dist[1]> 0.0008):
    if(temp_dist[0]< 0.0002 and temp_dist[1]> 0.0016):
        #print ('bbbbbbbbb')
        leastDistIndex = np.where(distArr == temp_dist[0])
        df_source['star_bool'][matchIndexArr[leastDistIndex[0][0]]] = 1
        
df_source.to_pickle(outFile)       
    




