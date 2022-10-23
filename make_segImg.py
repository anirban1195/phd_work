#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 10:49:03 2020

@author: anirban
"""

import numpy as np
from astropy.io import fits
n= 100
startArr = np.arange(0,100, 5)
meanArr  = np.array([100,200, 300, 400,500, 600])
output = '/scratch/halstead/d/dutta26/test_img/'
f=open('/scratch/halstead/d/dutta26/test_img/file_list.ascii', 'w+')
for j in range(n):
    print (j)
    a=np.zeros((5000, 5000), dtype =np.float32)
    a[:] = np.nan
    start = np.random.choice(startArr, 1)[0]
    start1 = np.random.choice(startArr, 1)[0]
    xPos = start
    yPos = start1
    mean = np.random.choice(meanArr, 1)[0]
    a[0:5000, 0:5000] = np.random.normal(mean, np.sqrt(mean), (5000,5000))
# =============================================================================
#     while(xPos< 5000 and yPos<5000):
#         a[yPos:yPos+1280, xPos:xPos+1280] = np.random.normal(mean, np.sqrt(mean), (1280,1280))
#         yPos = yPos + 1300
#         if(yPos > 5000):
#             yPos = start1
#             xPos = xPos + 1300
# =============================================================================
    
    hdu = fits.PrimaryHDU(a)  
    hdu.writeto(output + 'img_'+str(j)+'.fits', overwrite=True)
    f.write(output + 'img_'+str(j)+'.fits \n' )
    
f.close()
    
        
        