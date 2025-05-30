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
from scipy.fftpack import fftshift, ifftshift, fft2, ifft2

## polar coordinates
def Zernike_polar(coefficients, r, u):
    Z = coefficients
    Z1  =  Z[0] #piston
    Z2  =  Z[1]  * 2*r*np.cos(u) #tilt Y
    Z3  =  Z[2]  * 2*r*np.sin(u) #tilt X
    Z4  =  Z[3]  * np.sqrt(3)*(2*r**2-1) #defocus
    Z5  =  Z[4]  * np.sqrt(6)*r**2*np.sin(2*u) #astigX
    Z6  =  Z[5]  * np.sqrt(6)*r**2*np.cos(2*u) #astigY
    Z7  =  Z[6]  * np.sqrt(8)*(3*r**2-2)*r*np.sin(u) #coma X
    Z8  =  Z[7]  * np.sqrt(8)*(3*r**2-2)*r*np.cos(u) #coma Y
    
    ZW = Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8 
    return ZW

r_range= np.linspace(0,1,100)
theta_range= np.linspace(-np.pi, np.pi, 100)
coefficients=[0,0,0,0,10,0,0,0]

R, THETA= np.meshgrid(r_range, theta_range)

zernike= Zernike_polar(coefficients, R, THETA)

P= np.exp(1j*zernike) #Pupil function
PSF = ifftshift(fft2(fftshift(P)))
PSF= (np.abs(PSF))**2
X, Y = R*np.cos(THETA), R*np.sin(THETA)
#plt.contourf(X, Y, P, levels=100)
#plt.colorbar()
plt.imshow(PSF)
plt.show()
