#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:08:48 2024

@author: alvaro
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import projections
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import cumtrapz
from matplotlib.animation import FuncAnimation
import math
import sys


#Definition of functions and constants


#H = 10**(-1)                       #Hubble constant

Lzero= 0.5 
cf = 1/(4*np.pi*Lzero**2)               #Overall constant

#sq3 = math.sqrt(3)/(2*H)                  
N= int(sys.argv[1])
L = int(sys.argv[2])
R = int(sys.argv[3])
n_ea = int(sys.argv[4])
n_dh = int(sys.argv[5])

if len(sys.argv) != 6:
    sys.exit('Missing arguments, N, L, R, n_ea, n_dh required')

def switch(z, DH):               #Switching funcion
        return np.cos(np.pi*z/2*DH)**4 

def k(n, l, r):                           #k wavenumber
    return math.sqrt(n**2 + l**2 + r**2)

def nonunit(z, n, l, r):                     #Part of the standard coupling
    return np.exp(-z)/math.sqrt(k(n, l, r))

def RealDyn(z,EA, n, l, r):                   #Real part of the function to be integrated
    return np.cos(z*EA) * np.cos(2*np.pi * k(n, l, r) * np.exp(-z)/Lzero) + np.sin(z*EA) * np.sin(2*np.pi * k(n, l, r) * np.exp(-z)/Lzero) 

def ImDyn(z,EA, n, l, r):                     #Imaginary part of the function to be integrated
    return  np.cos(z*EA) * np.sin(2*np.pi * k(n, l, r) * np.exp(-z)/Lzero) - np.sin(z*EA) * np.cos(2*np.pi * k(n, l, r) * np.exp(-z)/Lzero) 

def fReal(z, DH, EA, n, l, r):                #We split the integration into the Real and Imaginary part
    return switch(z,DH) * nonunit(z, n, l, r) * RealDyn(z, EA, n, l, r)

def fImag(z, DH, EA, n, l, r):
    return switch(z,DH) * nonunit(z, n, l, r) * ImDyn(z, EA, n, l, r)


#Try graphic for specific values of EA and DH
#List with values of the function
#g_list= f(z_list,1,1)   #Arbitrary values for EA and DH
#plt.plot(z_list,g_list)

#Define Grid of parameters
EA_array = np.linspace(-12,12,n_ea)     #Space of parameters

print('len EA_array', len(EA_array))

DH_array = np.linspace(0.175,0.8,n_dh)    #Is good to take different lenghts so avoid confusion

print('len DH_array', len(DH_array))

X, Y = np.meshgrid(DH_array,EA_array) #Grid definition of parameters


#Canvas
fR_array=np.zeros( (len(EA_array),len(DH_array)) ) #Define a matrix of zeros = CANVAS (Rows, Columns)
fI_array=np.zeros( (len(EA_array),len(DH_array)) )
ResultadoIntegral = np.zeros( (len(EA_array),len(DH_array)), dtype=np.complex128 ) #We define a matrix of 
                                                        # complex zeros where we gonna join the solution
F_array=np.zeros( (len(EA_array),len(DH_array)) )       # We gonna use this matrix to turn everything real

#def calculate_integral(real, imaginary, ea, dh):
    
    
    
    #fR_array[ind_ea][ind_dh] + 1j * fI_array[ind_ea][ind_dh]
    

for ind_ea,ea in enumerate (EA_array):                  #For each definite value of EA
    for ind_dh,dh in enumerate (DH_array):              #For each definite value of DH
        z_list= np.linspace(-dh, dh, 100000)    #We define the limit of integration differently,
                                        # the integral is zero outside DH because of the switching function
        for n in range(0,N+1):
            for l in range(0, L+1):
                for r in range(0, R+1):
                    if n == l == r == 0:
                        continue
                    fR_array[ind_ea][ind_dh] = cumtrapz(fReal(z_list,ea,dh, n, l, r),z_list)[-1] #integration real part
                    fI_array[ind_ea][ind_dh] = cumtrapz(fImag(z_list,ea,dh, n, l, r),z_list)[-1] #integration imaginary part
                    ResultadoIntegral[ind_ea][ind_dh] += fR_array[ind_ea][ind_dh] + 1j * fI_array[ind_ea][ind_dh] #junction of the result
                    F_array[ind_ea][ind_dh] += np.abs(ResultadoIntegral[ind_ea][ind_dh])**2 * cf #We take the norm times the constant
                                                                        #Sum each integral until N
    print(ind_ea)                                       #We set a progress control

"""    
#Plot
#Multicolour
fig = plt.figure()                              #Define variable of plot
ax = fig.add_subplot(111, projection='3d')  
    #The figure is gonna be 3D

surface = ax.plot_surface(X, Y, F_array , cmap='viridis') #Colour plot
ax.view_init(elev=30, azim=45)
ax.set_xlabel(r'$H Delta$')                    #Labels 
ax.set_ylabel(r'$E H^{-1}$')                                              #r is for Latex
ax.set_zlabel('$F$')
plt.savefig(f"Lo={Lzero}, N={N}, L={L}, R={R}.png", dpi=650)
plt.show()
"""




#Monocolour
#ax = plt.axes(projection="3d")
#ax.plot_surface(X, Y, f_array)
#ax.set_xlabel(r'$\frac{E}{H}$')
#ax.set_ylabel(r'$\Delta H$')
#ax.set_zlabel('F')
#plt.show()


