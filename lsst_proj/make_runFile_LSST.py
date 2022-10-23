#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 16:52:30 2022

@author: dutta26
"""


import numpy as np
f=open('/home/dutta26/codes/lsst_proj/phosim_lsst2.sub', 'w+')
f.write('#!/bin/sh -l \n')
f.write('# FILENAME:  phosim_lsst2.sub \n')
f.write('cd /home/dutta26/apps/phosim-phosim_release-2479c7396c0c/  \n')
#f.write('cd /home/dutta26/Downloads/phosim_core/  \n')
for j in np.arange(70, 120):
    f.write('mkdir /scratch/halstead/d/dutta26/lsst/several_img_f2/100'+str(j)+'/  \n')
    f.write('./phosim /scratch/halstead/d/dutta26/lsst/lsst_catalog1/'+str(j)+'.txt_sheared -i lsst -c examples/quickbackground -o /scratch/halstead/d/dutta26/lsst/several_img_f2/100'
            +str(j)+'/ -w /scratch/halstead/d/dutta26/lsst/work1 \n')
    
    #f.write('mkdir /scratch/halstead/d/dutta26/lsst/3/100'+str(j)+'_flat/  \n')
    #f.write('./phosim /scratch/halstead/d/dutta26/lsst/wiyn_cat1/'+str(j)+'.txt3_flat -i wiyn_odi -c examples/quickbackground -o /scratch/halstead/d/dutta26/lsst/3/100'
    #        +str(j)+'_flat/ -w /scratch/halstead/d/dutta26/lsst/work3 \n')
f.close()