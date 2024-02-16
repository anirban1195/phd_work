#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:23:42 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np 
import os,sys
import subprocess
import gzip
import matplotlib.pyplot as plt
import shutil
from astropy.stats import sigma_clipped_stats


sexLoc = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/'
ref_catalog = '/scratch/bell/dutta26/wiyn_sim1/wcs_reference.cat'

def correct_wcs(idNo, filt, loc_img):
    #Open and read reference wcs catalog
    f=open(ref_catalog)
    cat = f.readlines()
    f.close()
    ref_ra =[]
    ref_dec = []
    ref_flux=[]
    count = 0
    
    #Read the list of ra dec and dec into arrays
    for j in range(len(cat)):
        if((cat[j].split()[0]) == '#'):
         continue
        
        if(float(cat[j].split()[1]) < 10000):
            continue
        ref_ra.append(float(cat[j].split()[5])) 
        ref_dec.append(float(cat[j].split()[6]))
        ref_flux.append(float(cat[j].split()[1]))
        count += 1
    #Keep the 200 brightest souce
    ref_ra = np.array(ref_ra)
    ref_dec = np.array(ref_dec)
    ref_flux = np.array(ref_flux)
    
    #Take the 200 brightest sources as reference (assumes all are at flux > 10k)
    if(len(ref_flux)> 200):
        topIndices = np.argpartition(ref_flux, -200)[-200:]
    else:
        print ('Unable to find enough stars')
        
    ref_ra = ref_ra[topIndices]
    ref_dec = ref_dec[topIndices]    
    
    
    #Reun sectractor on the folder
    file =loc_img+filt+'/'+str(idNo)+'/output.fits'    
    bashCommand = 'sex '+ file +' -DETECT_MINAREA 9 -DETECT_THRESH 5 -ANALYSIS_THRESH 5 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+loc_img+filt+'/'+str(idNo)+'/output.weight.fits'+' -CHECKIMAGE_TYPE NONE -CATALOG_NAME /scratch/bell/dutta26/wiyn_sim1/temp.cat'
    #bashCommand = 'sex '+ file +' -DETECT_MINAREA 5 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+loc_img+filt+'/'+str(idNo)+'/output.weight.fits'+' -CHECKIMAGE_TYPE NONE -CATALOG_NAME /scratch/bell/dutta26/wiyn_sim1/temp.cat'

    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
    output, error = process.communicate()
    
    #Read the seatractor output
    f=open('/scratch/bell/dutta26/wiyn_sim1/temp.cat')
    cat = f.readlines()
    f.close()
    raArr =[]
    decArr=[]
    fluxArr=[]
    
    #Read the ourput into array or ra and dec
    for j in range(len(cat)):
        if((cat[j].split()[0]) == '#'):
         continue
        
        if(float(cat[j].split()[1]) < 5000):
            continue
        raArr.append(float(cat[j].split()[5])) 
        decArr.append(float(cat[j].split()[6]))
        fluxArr.append(float(cat[j].split()[1]))
        
    fluxArr=np.array(fluxArr)
    raArr = np.array(raArr)
    decArr = np.array(decArr)    
    
    #Take the 100 brightest source
    if(len(fluxArr)> 100):
        topIndices = np.argpartition(fluxArr, -100)[-100:]
    else:
        print ('Unable to find stars')
        
    raArr = raArr[topIndices]
    decArr = decArr[topIndices]
    
    #Attempt to match ra and dec
    raShift=[]
    decShift=[]
    distShift=[]
    #Now match indices 
    for j in range(len(raArr)):
        ra = raArr[j]
        dec = decArr[j]
        delRa = (ref_ra - ra)
        delDec = (ref_dec - dec)
        dist = np.sqrt(delRa**2 + delDec**2)
        loc = np.where(dist <0.002000 )[0]
        #print (len(loc), loc)
        for l in range(len(loc)):
            raShift.append(delRa[loc][l])
            decShift.append(delDec[loc][l])
            distShift.append(dist[loc][l])
    
    
    #Find average shifts        
    avgRaShift =  sigma_clipped_stats(raShift)[1] 
    avgDecShift = sigma_clipped_stats(decShift)[1] 
    print (avgRaShift, avgDecShift)
    
    #Set the ra dec for main file
   
    f=fits.open(loc_img+filt+'/'+str(idNo)+'/output.fits' )
    hdr = (f[0].header)
    f.close()
    fits.setval(loc_img+filt+'/'+str(idNo)+'/output.fits' , 'CRVAL1', value=hdr['CRVAL1']+avgDecShift)
    fits.setval(loc_img+filt+'/'+str(idNo)+'/output.fits' , 'CRVAL2', value=hdr['CRVAL2']+avgRaShift)   
    #Set the ra dec for weight file
    f=fits.open(loc_img+filt+'/'+str(idNo)+'/output.weight.fits' )
    hdr = (f[0].header)
    f.close()
    fits.setval(loc_img+filt+'/'+str(idNo)+'/output.weight.fits' , 'CRVAL1', value=hdr['CRVAL1']+avgDecShift)
    fits.setval(loc_img+filt+'/'+str(idNo)+'/output.weight.fits' , 'CRVAL2', value=hdr['CRVAL2']+avgRaShift)   
    
    #Set ra dec for star file 
    f=fits.open(loc_img+filt+'/'+str(idNo)+"_star/output.fits")
    hdr = (f[0].header)
    f.close()
    fits.setval(loc_img+filt+'/'+str(idNo)+"_star/output.fits", 'CRVAL1', value=hdr['CRVAL1']+avgDecShift)
    fits.setval(loc_img+filt+'/'+str(idNo)+"_star/output.fits", 'CRVAL2', value=hdr['CRVAL2']+avgRaShift)   
    #Set ra dec for star weight file 
    f=fits.open(loc_img+filt+'/'+str(idNo)+"_star/output.weight.fits")
    hdr = (f[0].header)
    f.close()
    fits.setval(loc_img+filt+'/'+str(idNo)+"_star/output.weight.fits", 'CRVAL1', value=hdr['CRVAL1']+avgDecShift)
    fits.setval(loc_img+filt+'/'+str(idNo)+"_star/output.weight.fits", 'CRVAL2', value=hdr['CRVAL2']+avgRaShift)   
    
