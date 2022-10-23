#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 09:13:54 2022

@author: dutta26
"""


from astroquery.ned import Ned
from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy import coordinates
import numpy as np
Ned.TIMEOUT=500
co = coordinates.SkyCoord(ra=328.2833, dec=17.6669,
                          unit=(u.deg, u.deg), frame='icrs')
result_table = Ned.query_region(co, radius=0.9 * u.deg, equinox='J2000.0')
df = result_table.to_pandas()
df = df[df['Redshift'] >= 0]
raList = np.array(df['RA'])
decList = np.array(df['DEC'])


f=fits.open('/scratch/halstead/d/dutta26/abell_2390/redshift_map.fits', mode = 'update')
w = wcs.WCS(f[0].header)
data = np.array(f[0].data)
#First convert ra dec to x and y
tempList = np.zeros((len(raList), 2), dtype = np.float32)
tempList[:,0] = raList
tempList[:,1] = decList
temp = w.wcs_world2pix(tempList, 0)
xList = temp[:,0]
yList = temp[:,1]
xMax,yMax = np.shape(data)

for j in range(len(xList)):
    x=int(xList[j])
    y=int(yList[j])
    if(x>xMax or x<0 or y>yMax or y<0):
        continue
    f[0].data[y-25:y+25, x+25] = 100
    f[0].data[y-25:y+25,x-25] = 100
    f[0].data[y-25, x-25:x+25] = 100
    f[0].data[y+25, x-25:x+25] = 100
    
f.flush()
    
    
    
    
    
    
    
    
    
    
    
    
    
    