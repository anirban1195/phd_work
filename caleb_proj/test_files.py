#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:57:31 2020

@author: anirban
"""

from astropy.io import fits
import numpy as np
import os
import subprocess
import gzip
from astropy import units as u
from astropy.coordinates import SkyCoord
import time


phosimLoc = '/home/dutta26/apps/phosim-phosim_release-60db50aae970/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/halstead/d/dutta26/abell_2390/'
folder = ['g', 'r']
outFile= open('/home/dutta26/comparison_gr1.txt', 'w+')


for sets in os.listdir(dataSet):
    if(sets not in folder):
        continue
    for files in os.listdir(dataSet+sets):
        print (files)
        if('weight' in files):
            continue
