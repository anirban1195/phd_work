#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 07:27:17 2021

@author: anirban
"""

from astropy.io import fits
import numpy as np
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
#from astroquery.gaia import Gaia
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.coordinates import FK5 

raIndex = 5
decIndex =6
sizeIndex1 = 7
sizeIndex2 = 8
#Taken in a list of indices and returns the ones that are not very large 
def removeLarge(indices, catalog , psfSize = 3):
    with open(catalog) as f:
        content = f.readlines()
    returnIndexList =[]    
    for j in range(len(content)):
        flag = 0
        if ( j not in indices):
            continue
        if((content[j].split())[0] == '#'):
            continue
        temp = content[j].split()
        
        size1 = float(temp[sizeIndex1])
        size2 = float(temp[sizeIndex2])
        elongation = 1- min(size1, size2)/max(size1, size2)
        size = np.sqrt(size1*size2)
        if(size > 3*psfSize or elongation>0.9):
            flag = 1
            
        if(flag == 0):
            returnIndexList.append(j)
            
    return returnIndexList


def removeStars(indices, catalog):
    Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source" # Select early Data Release 3
    
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/abell_ir_coadd_wted1.fits')
    
    img = np.array(f[0].data)
    sizey, sizex = np.shape(img)
    w=wcs.WCS(f[0].header)
    f.close()
    del img
    [[startRa, startDec]]  = w.wcs_pix2world([[0,0]], 0)
    [[endRa, endDec]] = w.wcs_pix2world([[sizex, sizey]], 0)
    midRa = 0.5*(startRa +endRa )
    midDec = 0.5* (startDec + endDec)
    
    coord = SkyCoord(ra=midRa , dec=midDec, unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(np.abs(startRa - endRa ), u.deg)
    height = u.Quantity(np.abs(startDec - endDec), u.deg)
    Gaia.ROW_LIMIT = -1
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    r.pprint()
    
    df = r.to_pandas()
    raArr =[]
    decArr =[]
    indexArr =[]
    #Read the catalog 
    with open(catalog) as f:
        content = f.readlines()

    for j in range(len(content)):
        if((content[j].split())[0] == '#'):
           continue
        if(j not in indices):
            continue
        raArr.append(float(  (content[j].split())[raIndex]) )
        decArr.append(float(  (content[j].split())[decIndex]) )
        indexArr.append(j)
    
    
    returnIndexList =[]
    star_raArr = np.array(df['ra'])
    star_decArr = np.array(df['dec'])
    for j in range(len(raArr)):
        flag = 0
        print (j, len(df))
        ra = float( raArr[j])
        dec = float( decArr[j])
        for k in range(len(df)):
            ra_df = star_raArr[k]
            dec_df = star_decArr[k]
            
            if (abs(ra_df - ra)< 0.0005 and abs(dec_df - dec)< 0.0005):
                flag = 1
                break
        if(flag == 0):
            returnIndexList.append(indexArr[j])
            
    return returnIndexList



            