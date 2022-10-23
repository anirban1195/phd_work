#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:51:09 2022

@author: dutta26
"""
import numpy as np
import sys
f=open('/scratch/halstead/d/dutta26/lsst/lsst_catalog/lsst.cat')
content = f.readlines()
f.close()


z_mass = 0.2
cent_ra =  237.0083  
cent_dec = -15.2308  


f_new = open('/scratch/halstead/d/dutta26/lsst/lsst_catalog/sheared.cat', 'w+')
# =============================================================================
# f_new = open('/home/dutta26/apps/phosim-phosim_release-2479c7396c0c/output/sheared.cat', 'w+')
# f_new.write('rightascension 328.4083 \n')
# f_new.write('declination 17.6697 \n')
# f_new.write('altitude 67.52902083333335\n')
# f_new.write('azimuth 235.70144583333334\n')
# f_new.write('moonalt -9.976905131961816\n')
# f_new.write('moonaz 247.4479727986705\n')
# f_new.write('moonphase 55.654687364728574\n')
# f_new.write('sunalt -51.9601859976301\n')
# f_new.write('sunaz 359.7246167905777\n')
# f_new.write('seed 1001\n')
# f_new.write('vistime 60.0 \n')
# f_new.write('filter 4 \n')
# f_new.write('seeing 1.24938576\n')
# f_new.write('obshistid 1001 \n')
# arr=[]
# =============================================================================
for j in range(len(content)):
    if(len(content[j].split()) == 15):
        f_new.write(content[j])
    elif(len(content[j].split()) == 30):
        z= float((content[j].split())[6])
        
        #If in foreground
        if(z< z_mass):
            f_new.write(content[j])
            continue
        
        ra = float((content[j].split())[2])
        dec = float((content[j].split())[3])
        
        delRa = -(ra - cent_ra)*3600
        delDec = (dec - cent_dec)*3600
        
        r = np.sqrt(delRa**2 +delDec**2)
        cos2phi = (delRa**2 - delDec**2)/r**2
        sin2phi = 2*delRa*delDec/r**2
        
        shear_signal = 0.3/(r/60)**2
        if(shear_signal > 0.6):
            shear_signal = 0.6
        gamma1 = -shear_signal*cos2phi
        gamma2 = shear_signal*sin2phi
        
        for k in range(len(content[j].split())):
            if(k == 7):
                f_new.write(str(gamma1)+' ')
                continue
            if(k == 8):
                f_new.write(str(gamma2)+' ')
                continue
            f_new.write((content[j].split())[k] +' ')
        f_new.write('\n')
        
f_new.close()
        