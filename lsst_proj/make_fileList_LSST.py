#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 17:30:10 2022

@author: dutta26
"""


import os,shutil
outLoc = '/scratch/halstead/d/dutta26/lsst/lsst1/'
inp_file = '/scratch/halstead/d/dutta26/lsst/filter5/'
for file in os.listdir(inp_file):
    for file1 in os.listdir(inp_file+file+'/'):
        if(len(os.listdir(inp_file+file+'/')) == 153):
            if('C' in file1):
                continue
            else:
                print (inp_file+file+'/'+file1, outLoc+file1)
                shutil.copy(inp_file+file+'/'+file1, outLoc+file1)