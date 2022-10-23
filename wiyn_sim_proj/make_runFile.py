#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 09:41:20 2022

@author: dutta26
"""

import numpy as np
f=open('/home/dutta26/codes/wiyn_sim_proj/phosim_run3.sub', 'w+')
f.write('#!/bin/sh -l \n')
f.write('# FILENAME:  phosim_run3.sub \n')
f.write('cd /home/dutta26/Downloads/phosim_core/ \n')
for j in np.arange(65, 99):
    f.write('mkdir /scratch/halstead/d/dutta26/abell_2390/odi_img_realsee/100'+str(j)+'/  \n')
    f.write('./phosim examples/odi_catalog/'+str(j)+'.txt -i wiyn_odi -c examples/quickbackground -o /scratch/halstead/d/dutta26/abell_2390/odi_img_realsee/100'
            +str(j)+'/ -w /scratch/halstead/d/dutta26/abell_2390/work  \n')
f.close()