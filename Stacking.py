#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:47:37 2018

@author: jishnu
"""

"""
Requirements:
-------------
    Linux packages:
        rawtran

    Python packages:
        scipy
        numpy
        astropy
        matplotlib
        astroalign


ref :
    http://manpages.ubuntu.com/manpages/bionic/man1/rawtran.1.html
    https://www.numpy.org
    https://www.scipy.org
    https://www.astropy.org
    http://toros-astro.github.io/astroalign
"""

#%%
import os
import glob
import numpy as np
import astroalign as aa
from astropy.io import fits
import matplotlib.pyplot as plt

#%%

#os.chdir('/your/path/to/folder/with/images')

input_files = glob.glob('*.NEF') #CR2 for canon RAW files

#%%
def convert2FITS(RAWfile):
    output = str(RAWfile).replace('NEF','fits')
    print('converting %s to fits format'%RAWfile)

    os.system('rawtran '+RAWfile+' -o '+output)

    return output

#%%
def imgshow(image): #pass in an mxn array to be displayed
    median,std = np.median(image),np.std(image)
    plt.figure(figsize=(9,9)) #size
    plt.imshow(image,vmin=median-1*std,vmax=median+5*std,cmap='gray')
    #use the above line for contrast stretching
    plt.gca().invert_yaxis()
    plt.grid(False)
#    plt.colorbar()
    plt.show()
    plt.close() #clean up after yourself!!

#%%
def get_img(file): #returns 2D array from fits file data
    with fits.open(file) as hdul:
            image = hdul[0].data    #defining image
            image = np.array(image) #converting image to an array
    return image #returns a 3D data cube with 3 layers each for R,G,B colors

#%%
def align(image,refimage): #magic!!
    try:
        p, (pos_image,pos_refimage) = aa.find_transform(image,refimage)
        return p #transformation parameters
    except RuntimeWarning:
        image += abs(np.amin(image)+10)
#        image = np.log10(image)
        p, (pos_image,pos_refimage) = aa.find_transform(image,refimage)
        return p #transformation parameters

#%%
def stack(files,refimage):
    
    y,x = refimage.shape[-2],refimage.shape[-1]
    
    #creating empty arrays to hold each of RGB channels stacks
    stacksR = np.zeros((y, x, len(files)))
    stacksG = np.zeros((y, x, len(files)))
    stacksB = np.zeros((y, x, len(files)))
    
    rawHDU    = fits.open(files[0])
    rawHEADER = rawHDU[0].header
    
    for i in range(len(files)):
         
        rawDATA   = get_img(files[i]) #Image part of the file
        
        try:
            p          = align(rawDATA[0],refimage) #getting transformation
            rawDATAred = aa.apply_transform(p,rawDATA[0],refimage)
            rawDATAgrn = aa.apply_transform(p,rawDATA[1],refimage) 
            rawDATAblu = aa.apply_transform(p,rawDATA[2],refimage)
            
            stacksR[:,:,i] = rawDATAred
            stacksG[:,:,i] = rawDATAgrn
            stacksB[:,:,i] = rawDATAblu
        
        except:
            print('oops file %s didnt meet alignment criteria'%files[i])
            files.remove(files[i])
            continue
        
    reds = np.median(stacksR, axis=2)
    grns = np.median(stacksG, axis=2)
    blue = np.median(stacksB, axis=2)
    
    procDATA       = np.array([reds,grns,blue])
    procHDU        = fits.PrimaryHDU(procDATA) #new file, processed
    procHDU.header = rawHEADER #replacing header
    
    #specify your filename
    procHDU.writeto('stacked.fits', overwrite=True)
    
    
#%%

for file in input_files:
    output = convert2FITS(file)
    print('converted %s to fits'%output)
    
    
files = glob.glob('*.fits')

#%%
refimage = get_img(files[0]) # pick a reference image with respect to which 
                                   # will align all other images

stack(files,refimage[0])


#%%
'''
import rawpy
refimage = rawpy.imread(input_files[0])

from scipy.ndimage.filters import gaussian_filter

#%%
rgb = refimage.postprocess(gamma=(1,1),
                           no_auto_bright=True, 
                           output_bps=16)


imgshow(gaussian_filter(rgb[:,:,1][1500:3000,2500:4000], sigma=3))
#%%
plt.imshow(rgb[:,:,0])
plt.colorbar()


#%%
def get_img_fromraw(file):
    
    raw = rawpy.imread(file)
    rgb = raw.postprocess(gamma=(1,1),no_auto_bright=True, output_bps=16)
    return rgb
    
    
#%%
def stack_RAW(input_files,refimage):
    
    y,x = refimage.shape[-2],refimage.shape[-1]
    
    #creating empty arrays to hold each of RGB channels stacks
    stacksR   = np.zeros((y, x, len(input_files)))
    stacksG   = np.zeros((y, x, len(input_files)))
    stacksB   = np.zeros((y, x, len(input_files)))
    
#    raw       = rawpy.imread(input_files[0])
#    rgb       = raw.postprocess(gamma=(1,1),
#                                     no_auto_bright=True, 
#                                     output_bps=16)
    
    rawHDU    = fits.open(files[0])
    rawHEADER = rawHDU[0].header
    
    for i in range(len(input_files)):
         
        rawDATA   = get_img_fromraw(input_files[i]) #Image part of the file
        
        try:
            p          = align(rawDATA[:,:,1],refimage) #getting transformation
            rawDATAred = aa.apply_transform(p,rawDATA[:,:,0],refimage)
            rawDATAgrn = aa.apply_transform(p,rawDATA[:,:,1],refimage) 
            rawDATAblu = aa.apply_transform(p,rawDATA[:,:,2],refimage)
            
            stacksR[:,:,i] = rawDATAred
            stacksG[:,:,i] = rawDATAgrn
            stacksB[:,:,i] = rawDATAblu
        
        except:
            print('oops file %s didnt meet alignment criteria'%files[i])
            files.remove(files[i])
            continue
        
    reds = np.median(stacksR, axis=2)
    grns = np.median(stacksG, axis=2)
    blue = np.median(stacksB, axis=2)
    
    procDATA       = np.array([reds,grns,blue])
    procHDU        = fits.PrimaryHDU(procDATA) #new file, processed
    procHDU.header = rawHEADER #replacing header
    
    #specify your filename
    procHDU.writeto('stacked2.fits', overwrite=True)
    
#%%
#input_files = glob.glob('*.NEF') 

#ref = get_img_fromraw(input_files[0])[:,:,1]
#stack_RAW(input_files,ref)
'''

