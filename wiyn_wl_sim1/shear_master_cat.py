#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 12:09:02 2023

@author: dutta26
"""

import numpy as np
import sys
import getgamma 
f=open('/scratch/bell/dutta26/wiyn_sim1/master.cat')
content = f.readlines()
f.close()


z_mass = 0.0
cent_ra =  328.3941
cent_dec = 16.6697


f_new = open('/scratch/bell/dutta26/wiyn_sim1/sheared_uniform.cat', 'w+')
for j in range(len(content)):
    if(len(content[j].split()) == 15):
        f_new.write(content[j])
    elif(len(content[j].split()) == 30):
        zs= float((content[j].split())[6])
        
        #If in foreground
        if(zs< z_mass):
            f_new.write(content[j])
            continue
        
        ra = float((content[j].split())[2])
        dec = float((content[j].split())[3])
        rs = 0.25
        #gamma1, gamma2, kappa = getgamma.getGamma(cent_ra, cent_dec, z_mass, ra, dec, zs, rs)
        gamma1 = 0.0
        gamma2 =0.05
        kappa =0
        for k in range(len(content[j].split())):
            if(k == 7):
                f_new.write(str(gamma1)+' ')
                continue
            if(k == 8):
                f_new.write(str(gamma2)+' ')
                continue
            if(k == 9):
                f_new.write(str(kappa)+' ')
                continue
            f_new.write((content[j].split())[k] +' ')
        f_new.write('\n')
        
f_new.close()
        