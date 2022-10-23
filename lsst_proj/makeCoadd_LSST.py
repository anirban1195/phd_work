#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 08:03:28 2022

@author: dutta26
"""


import os,sys
from astropy.io import fits 
import subprocess
swarpLoc = '/home/dutta26/apps/bin/bin/'

#filt = 'filter1'
filtList= ['filter0', 'filter1', 'filter2', 'filter3','filter4','filter5']
for filt in filtList:
    loc = '/scratch/halstead/d/dutta26/lsst/'+filt+'/'
    f=open('/home/dutta26/lsst_temp.ascii', 'w+' )
    for folder in os.listdir(loc):
        if('coadd' in folder):
            continue
        if('flat' in folder):
            continue
        for sub_folder in os.listdir(loc+folder):
            if('coadd' not in sub_folder):
                continue
            if('weight' in sub_folder):
                continue
            f.write(loc+folder+'/'+sub_folder+'\n')
            print (sub_folder)
           
    f.close()     
    #sys.exit()
    bashCommand = './swarp @/home/dutta26/lsst_temp.ascii -c /home/dutta26/default1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/lsst/'+filt+'/coadd_all.fits -WEIGHTOUT_NAME /scratch/halstead/d/dutta26/lsst/'+filt+'/coadd_all.weight.fits -RESAMPLE_DIR /scratch/halstead/d/dutta26/lsst/ '
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
    
    f=open('/home/dutta26/lsst_temp.ascii', 'w+' )
    for folder in os.listdir(loc):
        if('coadd' in folder):
            continue
        if('flat' not in folder):
            continue
        for sub_folder in os.listdir(loc+folder):
            if('coadd' not in sub_folder):
                continue
            if('weight' in sub_folder):
                continue
            f.write(loc+folder+'/'+sub_folder+'\n')
            print (sub_folder)
           
    f.close()   
    #sys.exit()  
    bashCommand = './swarp @/home/dutta26/lsst_temp.ascii -c /home/dutta26/default1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/lsst/'+filt+'/coadd_flat_all.fits -WEIGHTOUT_NAME /scratch/halstead/d/dutta26/lsst/'+filt+'/coadd_flat_all.weight.fits -RESAMPLE_DIR /scratch/halstead/d/dutta26/lsst/ '
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
