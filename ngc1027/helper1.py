#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:08:46 2020

@author: dutta26
"""


import numpy as np
from scipy.optimize import curve_fit

def genPoly( xy_tuple, A, B,C,D,E,F, G, H):
    (x , y) = xy_tuple
    a= A*x**3 + B*y**3 + C*x**2 + D*y**2 + E*x*y**2 + F*y*x**2 + G*x*y + H
    #print (np.shape(a), np.shape(x))
    return a.ravel()

def returnPolyVal(popt, x, y):
    return popt[0]*x**3 + popt[1]*y**3 + popt[2]*x**2 + popt[3]*y**2 + popt[4]*x*y**2 + popt[5]*y*x**2 + popt[6]*x*y + popt[7]   

def fitPoly2D(X,Y, Z, sigma):
    # Initial guesses to the fit parameters.
    p0 = (0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
                  
    
    # Do the fit, using our custom _gaussian function which understands our
    # flattened (ravelled) ordering of the data points.
    popt, pcov = curve_fit(genPoly, (X,Y) , Z, sigma=sigma)
    return popt,pcov