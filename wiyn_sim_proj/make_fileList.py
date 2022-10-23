#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 09:21:51 2022

@author: dutta26
"""


import os,shutil
outLoc = '/scratch/halstead/d/dutta26/abell_2390/odi_img1/'
inp_file = '/scratch/halstead/d/dutta26/abell_2390/np_data/'
count = 0
for file in os.listdir(inp_file):
    print (len(os.listdir(inp_file+file+'/')), file)
    if(len(os.listdir(inp_file+file+'/')) == 1690):
        count += 1
        for file1 in os.listdir(inp_file+file+'/'):
            if('CELL' in file1 ):
                continue
            else:
                #print (inp_file+file+'/'+file1, outLoc+file1)
                shutil.copy(inp_file+file+'/'+file1, outLoc+file1)
    if(count >= 200):
        break
            
    
