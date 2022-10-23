#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 19:08:23 2022

@author: dutta26
"""

store = np.load()

cat = f.open()
content = f.readlines()
f.close()

store_cat = np.zeros((len(content),5), dtype = np.float32 )
count = 0
for j in range(len(content)):
    
    z= float((content[j].split())[6])
    ra = float((content[j].split())[2])
    dec = float((content[j].split())[3])
    mag = float((content[j].split())[4])
    obj = (content[j].split())[12] 
    if(mag>29):
        continue
    id_num = (content[j].split()[1]).split('.')[0]
    #store_cat[j,0:4] = ra,dec,mag,z
    if('sersic' in obj):
        prev_id_num = (content[j-1].split()[1]).split('.')[0]
    if(id_num == prev_id_num):
        continue
    
    store_cat[j,0:4] = ra,dec,mag,z
    count += 1