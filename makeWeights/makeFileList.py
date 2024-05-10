import os
import sys
from astropy.io import fits
loc = str(sys.argv[1])
#loc = '/scratch/halstead/d/dutta26/abell_2390/i/'
#lensQuality = int(sys.argv[3])
lensQuality = 0
outLoc = str(sys.argv[2])
#outLoc = '/home/dutta26/codes/makeWeights/fileList_swarp_i.ascii'
f=open(outLoc, 'w+')
#print (sys.argv[1])
for files in os.listdir(loc):
#This is for 2390
    print (files)
    if('temp' in files):
        continue
    f1=fits.open(loc+files)
    back = float((f1[0].header)['SKYBG'])
    seeing = float((f1[0].header)['SEEING'])
    f1.close()
    
    if(lensQuality == 1 and seeing>1.2): #then do not use the image
        continue
    
    elif('.weight' in files):
        continue
    
    else:
        #This is for 2390
        if('/r/' in loc and ('84728.2' in files or '90144.8' in files or '03729.4' in files or
                             '92055.5' in files or '92055.6' in files or '92055.9' in files or
                             '4940.1' in files)):
            continue
        if('/i/' in loc and ('23106.2' in files or '01046.5' in files or '83209' in files or
                             '20044' in files or '0611.7' in files or '0611.6' in files or '0611.3' in files
                             or '0136.2' in files or '15637.6' in files) ):
            continue
        
        if('/g/' in loc and ('1110.2' in files or '5530.1' in files or '4440.9' in files or '3024.5' in files or
                             '2605.6' in files or '903.4' in files or '4903.3' in files or '4903.2' in files or 
                             '4903.1' in files or '5951.1' in files or '1346.1' in files or '1018.1' in files or 
                             '0507.1' in files or '5756.1' in files or '2642.1' in files or '84517.7' in files)):
            continue
        
        if('/z/' in loc and ('1318.6' in files or '2722.2' in files or '3256.8' in files )):
            continue
        
        if('/u/' in loc and ('2044.9' in files )):
            continue
        
        f.write(loc+files+'\n')
# =============================================================================
#     #This is for 2219
#     if('.weight' in files):
#         continue
#     else:
#         if('/r/' in loc and ('4824.1' in files or '2925.6' in files or '2925.3' in files or
#                              '5734.4' in files)):
#             continue
#         if('/i/' in loc and ('11226.6' in files )):
#             continue
#         if('/g/' in loc and ('02630.1' in files )):
#             continue
#         
#         f.write(loc+files+'\n')
# =============================================================================
		
    
f.close()
