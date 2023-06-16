#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 06:52:02 2023

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt
a = np.load('/home/dutta26/codes/sersic_test/data7.npy')
fmtArr = ['r.', 'b.', 'k.', 'g.']

for j in range(len(a)):
    plt.errorbar(a[2,j], a[0,j], yerr =a[1,j] , xerr = a[3,j], fmt = fmtArr[j])