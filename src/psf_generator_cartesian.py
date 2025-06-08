#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Wed May 28 11:59:16 PM IST 2025
@author: janmejoy
@hostname: machine

DESCRIPTION
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2

## polar coordinates
def Zernike_polar(coefficients, x, y):
    r, u= (x**2+y**2), np.atan2(y,x)
    Z = coefficients
    Z0  =  Z[0] #piston
    Z1  =  Z[1]  * 2*r*np.cos(u) #tilt X
    Z2  =  Z[2]  * 2*r*np.sin(u) #tilt Y
    Z3  =  Z[3]  * np.sqrt(3)*(2*r**2-1) #defocus
    Z4  =  Z[4]  * np.sqrt(6)*r**2*np.cos(2*u) #astigX
    Z5  =  Z[5]  * np.sqrt(6)*r**2*np.sin(2*u) #astigY
    Z6  =  Z[6]  * np.sqrt(8)*(3*r**2-2)*r*np.cos(u) #coma X
    Z7  =  Z[7]  * np.sqrt(8)*(3*r**2-2)*r*np.sin(u) #coma Y
    Z8  =  Z[8]  * np.sqrt(5)*(6*r**4-6*r**2+1) #primaryspherical
    Z9  =  Z[9]  * np.sqrt(8)*r**3*np.cos(3*u) #trefoilX
    Z10 =  Z[10] * np.sqrt(8)*r**3*np.sin(3*u) #trefoilY
    
    ZW = Z0+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10 
    return ZW

def psf(coefficients):
    _x=_y= np.linspace(-1,1,N)
    X,Y= np.meshgrid(_x,_y)
    zernike= Zernike_polar(coefficients, X, Y)
    mask = (X**2+Y**2)< 1 # make a circular aperture of unit radius
    zernike= zernike*mask
    P= np.exp(-1j*zernike) #Complex Pupil Function
    PSF = ifftshift(fft2(fftshift(P)))
    #fftshift shifts the origin to center for computing FFT.
    #ifftshift shifts the origin back to the corner after performing FFT.
    PSF= (np.abs(PSF))**2
    PSF= PSF/np.sum(PSF) # Power normalized PSF
    return (zernike, PSF)
    
def visualize(zernike, PSF):
    ## Visualizations ##
    extent_limits= np.array([-frame_size, frame_size, -frame_size, frame_size]) #arcsec 
    plt.figure()
    plt.subplot(1,2,1)
    plt.title('Wavefront Map')
    plt.imshow(zernike, cmap='jet', origin='lower')
    plt.subplot(1,2,2)
    plt.title('PSF')
    plt.imshow(PSF, extent=extent_limits, cmap='RdBu', origin='lower', norm=LogNorm())
    plt.show()

if __name__=='__main__':
    N=4096
    D= 150e3 #u m
    frame_size= D*0.7/(12) #arcsec
    coefficients=[0 ,#piston 
                  0 ,#tiltX
                  0 ,#tiltY
                  -19.368*2 ,#defocus
                  -18.730*2 ,#astigX
                  -15.694*2 ,#astigY
                  -26.257*2 ,#comaX
                  -10.104*2 ,#comaY
                  -13.916*2 ,#primaryspherical
                  -20.271*2 ,#trefoilX
                  -18.015*2  #trefoilY 
                  ]
    gen_zernike, gen_PSF= psf(coefficients)
    visualize(gen_zernike, gen_PSF)
    
    '''
    _i= _j= np.arange(N)
    I,J= np.meshgrid(_i, _j)
    circ= np.sqrt((I-N/2)**2+(J-N/2)**2)
    ee_list=[]
    for rad in np.arange(100):
        circ_mask= circ<rad
        ee= np.sum(gen_PSF*circ_mask)
        ee_list.append(ee)

    plt.figure("EE plot")
    plt.plot(ee_list)
    plt.show()
    '''
