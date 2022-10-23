#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 22:14:23 2022

@author: dutta26
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys,gzip,shutil
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess
import zipfile
import glob

catalogLoc = '/scratch/halstead/d/dutta26/lsst/lsst_catalog1/0.txt_sheared'
#Read catalog 
f=open(catalogLoc)
content = f.readlines()
f.close()
star_ra =[]
star_dec = []
star_mag=[]
count = 0

for j in range(len(content)):
    if(len(content[j].split()) == 15):
        if(float(content[j].split()[4]) > 23): #filt 0 specific
            continue
        star_ra.append(   float(content[j].split()[2]) )
        star_dec.append(   float(content[j].split()[3]) )
        star_mag.append(   float(content[j].split()[4]) )
        count += 1

star_ra = np.array(star_ra)
star_dec = np.array(star_dec)

star_mag = np.array(star_mag)






# =============================================================================
# f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp_u.cat')
# content = f.readlines()
# f.close()
# star_ra =[]
# star_dec = []
# for j in range(len(content)):
#     if((content[j].split()[0]) == '#'):
#      continue
#     
#     star_ra.append(float(content[j].split()[5])) 
#     star_dec.append(float(content[j].split()[6]))
# =============================================================================




catalogLoc_flat = '/scratch/halstead/d/dutta26/lsst/lsst_catalog1/0.txt_flat'
#Read catalog 
f=open(catalogLoc_flat)
content = f.readlines()
f.close()
star_ra_flat =[]
star_dec_flat = []
star_mag_flat=[]
count = 0
for j in range(len(content)):
    if(len(content[j].split()) == 15):
        if(float(content[j].split()[4]) > 23):
            continue
        star_ra_flat.append(   float(content[j].split()[2]) )
        star_dec_flat.append(   float(content[j].split()[3]) )
        star_mag_flat.append(   float(content[j].split()[4]) )
        count += 1

star_ra_flat = np.array(star_ra_flat)
star_dec_flat = np.array(star_dec_flat)


star_thresh = 500 #filt 0 specific
folder ='/scratch/halstead/d/dutta26/lsst/'
outLoc= '/scratch/halstead/d/dutta26/lsst/lsst1/'
swarpLoc = '/home/dutta26/apps/bin/bin/'
sexLoc = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/'
destfolderLoc = '/scratch/halstead/d/dutta26/lsst/lsst2/'
for sub_folder in os.listdir(folder):
    if('filter' not in sub_folder):
        continue
    if('filter0' in sub_folder):
        star_thresh = 500
    else:
        star_thresh = 1000
    
    for sub_folder1 in os.listdir(folder+sub_folder):
        if('10' not in sub_folder1):
            continue
        
        #if('flat' not in sub_folder1):
        #    continue
         #Clear lsst and lsst2
        files = glob.glob('/scratch/halstead/d/dutta26/lsst/lsst1/*')
        for f in files:
            os.remove(f)
            
        files = glob.glob('/scratch/halstead/d/dutta26/lsst/lsst2/*')
        for f in files:
            os.remove(f)
                    
        
        print ('aaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
        
        
        for sub_folder2 in os.listdir(folder+sub_folder+'/'+sub_folder1):
            if('coadd' in sub_folder2 or 'sample' in sub_folder2):
                #os.remove(folder+sub_folder+'/'+sub_folder1+'/'+sub_folder2)
                continue
            if('C' in sub_folder2):
                continue
            else:
                print (folder+sub_folder+'/'+sub_folder1+'/'+sub_folder2, outLoc+sub_folder2)
                shutil.copy(folder+sub_folder+'/'+sub_folder1+'/'+sub_folder2, outLoc+sub_folder2)
                
        #Unzip the files 
        for sub_folder3 in os.listdir(outLoc):
            if('gz' not in sub_folder3):
                continue
            with gzip.open(outLoc+sub_folder3, 'rb') as f_in:
                with open(outLoc+sub_folder3[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
        
                    
        
        if('flat' in sub_folder1):
            
            #Unzip the files 
            for sub_folder3 in os.listdir(outLoc):
                if('gz' in sub_folder3):
                    continue
                bashCommand = './sex '+ outLoc+sub_folder3
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
                output, error = process.communicate()
                
                f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp.cat')
                cat = f.readlines()
                f.close()
                raArr =[]
                decArr=[]
                fluxArr=[]
                for j in range(len(cat)):
                    if((cat[j].split()[0]) == '#'):
                     continue
                    
                    if(float(cat[j].split()[1]) < star_thresh):
                        continue
                    raArr.append(float(cat[j].split()[5])) 
                    decArr.append(float(cat[j].split()[6]))
                    fluxArr.append(float(cat[j].split()[1]))
                
                fluxArr=np.array(fluxArr)
                raArr = np.array(raArr)
                decArr = np.array(decArr)    
                if(len(fluxArr)> 80):
                    topIndices = np.argpartition(fluxArr, -80)[-80:]
                elif(len(fluxArr) > 25):
                    topLen = len(fluxArr) - 2
                    topIndices = np.argpartition(fluxArr, -topLen)[-topLen:]
                else:
                    continue
                    
                raList = raArr[topIndices]
                decList = decArr[topIndices]
                #sys.exit()
                
                raShift=[]
                decShift=[]
                distShift=[]
                #Now match indices 
                for j in range(len(raList)):
                    ra = raList[j]
                    dec = decList[j]
                    delRa = (star_ra_flat - ra)
                    delDec = (star_dec_flat - dec)
                    dist = np.sqrt(delRa**2 + delDec**2)
                    loc = np.where(dist <0.001000 )[0]
                    print (len(loc), loc)
                    for l in range(len(loc)):
                        raShift.append(delRa[loc][l])
                        decShift.append(delDec[loc][l])
                        distShift.append(dist[loc][l])
                        
                if(len(raShift)< 20):
                    continue
                avgRaShift =  sigma_clipped_stats(raShift)[1] 
                avgDecShift = sigma_clipped_stats(decShift)[1] 
                #sys.exit()
                print (sigma_clipped_stats(raShift))
                print (sigma_clipped_stats(decShift))
                shutil.copy(outLoc+sub_folder3, destfolderLoc+sub_folder3)
                f=fits.open(destfolderLoc+sub_folder3)
                hdr = (f[0].header)
                f.close()
                
                
                print (hdr['CRVAL1'],avgRaShift, hdr['CRPIX1'])
                fits.setval(destfolderLoc+sub_folder3, 'CRVAL1', value=hdr['CRVAL1']+avgRaShift)
                fits.setval(destfolderLoc+sub_folder3, 'CRVAL2', value=hdr['CRVAL2']+avgDecShift)
                #sys.exit()
            
            #Now coadd
            f_coadd = open('/home/dutta26/lsst_temp.ascii', 'w+' )
            for sub_folder3 in os.listdir(destfolderLoc):
                f_coadd.write(destfolderLoc+ sub_folder3+'\n')
            f_coadd.close()
                
            
            bashCommand = './swarp @/home/dutta26/lsst_temp.ascii -c /home/dutta26/default1.swarp -IMAGEOUT_NAME '+folder+sub_folder+'/'+sub_folder1+'/final_coadd.fits -WEIGHTOUT_NAME '+folder+sub_folder+'/'+sub_folder1+'/final_coadd.weight.fits -RESAMPLE_DIR /scratch/halstead/d/dutta26/lsst/' 
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
            output, error = process.communicate()
        
        
            for sub_folder3 in os.listdir('/scratch/halstead/d/dutta26/lsst/lsst2/'):
                shutil.copy('/scratch/halstead/d/dutta26/lsst/lsst2/'+sub_folder3, folder+sub_folder+'/'+sub_folder1+'/sample.fits')
            
                
            
        else:
            
            #Unzip the files 
            for sub_folder3 in os.listdir(outLoc):
                if('gz' in sub_folder3):
                    continue
                bashCommand = './sex '+ outLoc+sub_folder3
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
                output, error = process.communicate()
                
                
                f=fits.open(outLoc+sub_folder3)
                hdr = f[0].header
                f.close()
                crpix1 = float(hdr['CRPIX1'])
                crpix2 = float(hdr['CRPIX1'])
                crval1 = float(hdr['CRVAL1'])
                crval2 = float(hdr['CRVAL2'])
                cd11 =float(hdr['CD1_1'])
                cd12 =float(hdr['CD1_2'])
                cd21 =float(hdr['CD2_1'])
                cd22 =float(hdr['CD2_2'])
                
                angle =  np.arctan2(cd21, cd11)
                pixScale = np.sqrt(cd11**2 +cd22**2) 
                
                f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp.cat')
                cat = f.readlines()
                f.close()
                raArr =[]
                decArr=[]
                fluxArr=[]
                for j in range(len(cat)):
                    if((cat[j].split()[0]) == '#'):
                     continue
                    
                    if(float(cat[j].split()[1]) < star_thresh):
                        continue
                    raArr.append(float(cat[j].split()[5])) 
                    decArr.append(float(cat[j].split()[6]))
                    fluxArr.append(float(cat[j].split()[1]))
                
                fluxArr=np.array(fluxArr)
                raArr = np.array(raArr)
                decArr = np.array(decArr)    
                if(len(fluxArr)> 80):
                    topIndices = np.argpartition(fluxArr, -80)[-80:]
                elif(len(fluxArr) > 55):
                    topLen = len(fluxArr) - 2
                    topIndices = np.argpartition(fluxArr, -topLen)[-topLen:]
                else:
                    continue
                #continue   
                raList = raArr[topIndices]
                decList = decArr[topIndices]
                #sys.exit()
                
                X = []
                Y =[]
                totMatched =0
                #Now match indices 
                for j in range(len(raList)):
                    ra = raList[j]
                    dec = decList[j]
                    delRa = (star_ra - ra)
                    delDec = (star_dec - dec)
                    dist = np.sqrt(delRa**2 + delDec**2)
                    loc = np.where(dist <0.001000 )[0]
                    print (len(loc), loc)
                    
                            
                    if(len(loc)> 0):
                        totMatched += 1
                        bLoc = loc[0]
                        bMag = star_mag[loc[0]]
                        for idx in loc:
                            if(bMag> star_mag[idx]):
                                bLoc = idx
                                bMag = star_mag[idx]
                        X.append([ ra - crval1, -dec+crval2, 1, 0])
                        X.append([ dec-crval2, ra - crval1, 0, 1])
                        Y.append([star_ra[bLoc]- crval1 ])
                        Y.append([star_dec[bLoc] - crval2])
                        
                        
                
                if(totMatched < 25):
                    continue
                #sys.exit()
                X = np.array(X)
                Y=np.array(Y)
                Z = np.linalg.lstsq(X,Y)[0]
                avgRaShift =  Z[2][0]
                avgDecShift = Z[3][0]
                avgAngleShift = np.arcsin(Z[1][0])
                #sys.exit()
                print (avgRaShift, avgDecShift, avgAngleShift)
                shutil.copy(outLoc+sub_folder3, destfolderLoc+sub_folder3)
                f=fits.open(destfolderLoc+sub_folder3)
                hdr = (f[0].header)
                f.close()
                
                
                print (hdr['CRVAL1'],avgRaShift, hdr['CRPIX1'])
                print (cd11, 5.55e-5*np.cos(angle+avgAngleShift))
                print (cd12 , -5.55e-5*np.sin(angle+avgAngleShift))
                print (cd21, 5.55e-5*np.sin(angle+avgAngleShift))
                print (cd22, 5.55e-5*np.cos(angle+avgAngleShift))
                #sys.exit()
                fits.setval(destfolderLoc+sub_folder3, 'CRVAL1', value=hdr['CRVAL1']+avgRaShift)
                fits.setval(destfolderLoc+sub_folder3, 'CRVAL2', value=hdr['CRVAL2']+avgDecShift)
                fits.setval(destfolderLoc+sub_folder3, 'CD1_1', value= 5.55e-5*np.cos(angle+avgAngleShift))
                fits.setval(destfolderLoc+sub_folder3, 'CD1_2', value= -5.55e-5*np.sin(angle+avgAngleShift))
                fits.setval(destfolderLoc+sub_folder3, 'CD2_1', value= 5.55e-5*np.sin(angle+avgAngleShift))
                fits.setval(destfolderLoc+sub_folder3, 'CD2_2', value= 5.55e-5*np.cos(angle+avgAngleShift))
                #sys.exit()
            
            #Now coadd
            f_coadd = open('/home/dutta26/lsst_temp.ascii', 'w+' )
            for sub_folder3 in os.listdir(destfolderLoc):
                f_coadd.write(destfolderLoc+ sub_folder3+'\n')
            f_coadd.close()
             
# =============================================================================
#             print (len(os.listdir(destfolderLoc)))
#             if(len(os.listdir(destfolderLoc)) < 9) :
#                 arr=['1008', '1009']
#                 if(sub_folder1 in arr ):
#                     continue
#                 sys.exit()
# =============================================================================
            bashCommand = './swarp @/home/dutta26/lsst_temp.ascii -c /home/dutta26/default1.swarp -IMAGEOUT_NAME '+folder+sub_folder+'/'+sub_folder1+'/final_coadd.fits -WEIGHTOUT_NAME '+folder+sub_folder+'/'+sub_folder1+'/final_coadd.weight.fits -RESAMPLE_DIR /scratch/halstead/d/dutta26/lsst/' 
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
            output, error = process.communicate()
        
        
            for sub_folder3 in os.listdir('/scratch/halstead/d/dutta26/lsst/lsst2/'):
                shutil.copy('/scratch/halstead/d/dutta26/lsst/lsst2/'+sub_folder3, folder+sub_folder+'/'+sub_folder1+'/sample.fits')
            
        #Clear lsst and lsst2
        files = glob.glob('/scratch/halstead/d/dutta26/lsst/lsst1/*')
        for f in files:
            os.remove(f)
            
        files = glob.glob('/scratch/halstead/d/dutta26/lsst/lsst2/*')
        for f in files:
            os.remove(f)
            
        #sys.exit()
            
            