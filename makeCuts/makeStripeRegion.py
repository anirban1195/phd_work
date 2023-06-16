#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 16:50:54 2022

@author: dutta26
"""
import numpy as np
temp =[]
for j in range(1,2):
    for k in range(2,3):
        for l in range(7):
            temp.append([j,k,l])


for j in range(5,6):
    for k in range(0,1):
        for l in range(7):
            temp.append([j,k,l])


for j in range(6,7):
    for k in range(1,2):
        for l in [3,5]:
            temp.append([j,k,l])
            

for j in range(8,9):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])
            
            
for j in range(9,10):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])
            
            
for j in range(11,12):
    for k in range(0,3):
        for l in range(7):
            temp.append([j,k,l])
            
            
for j in range(14,13):
    for k in range(0,3):
        for l in range(7):
            temp.append([j,k,l])            
            


for j in range(15,16):
    for k in range(0,3):
        for l in range(7):
            temp.append([j,k,l])


for j in range(16,17):
    for k in range(0,1):
        for l in range(7):
            temp.append([j,k,l])



for j in range(17,18):
    for k in [0,2,3,5]:
        for l in range(7):
            temp.append([j,k,l])
            

for j in range(21,22):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])



for j in range(23,24):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])


for j in range(24,25):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])


for j in range(25,26):
    for k in range(0,2):
        for l in range(7):
            temp.append([j,k,l])


stripeLocs = np.ones((31,8,8), dtype= np.float32)
for j in temp:
    print (j)
    stripeLocs[j[0], j[1], j[2]] =0.0

np.save('/home/dutta26/codes/stripeLoc.npy', stripeLocs)









