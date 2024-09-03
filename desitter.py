#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 18:08:48 2024

@author: alvaro
"""
import argparse
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
Lzero= 0.5 
cf = 1/(4*np.pi*Lzero**2)               #Overall constant

#-----------------------------------------------------------------------------
#Conventional Coupling
#-----------------------------------------------------------------------------

def switch(*, z_list, DH):               #Switching funcion
    return np.cos(np.pi * z_list/(2 * DH)) ** 4

def k(*, n, l, r):                           #k wavenumber
    return math.sqrt(n**2 + l**2 + r**2)

def overall(*, z_list, DH, n, l, r):                     #Part of the standard coupling
    return switch(z_list=z_list, DH=DH) * np.exp(-z_list)/(Lzero * math.sqrt(k(n=n, l=l ,r=r)))

def fReal1(*, z_list, EA, n, l, r):                #We split the integration into the Real and Imaginary part
    return np.cos(z_list * EA) * np.cos(2 * np.pi * np.exp(-z_list)/Lzero * k(n=n, l=l ,r=r))

def fReal2(*, z_list, EA, n, l, r):
    return np.sin(z_list * EA) * np.sin(2 * np.pi * np.exp(-z_list)/Lzero * k(n=n, l=l ,r=r))

def fImag1(*, z_list, EA, n, l, r):
    return np.cos(z_list * EA) * np.sin(2 * np.pi * np.exp(-z_list)/Lzero * k(n=n, l=l ,r=r))

def fImag2(*, z_list, EA, n, l, r):
    return np.sin(z_list * EA) * np.cos(2 * np.pi * np.exp(-z_list)/Lzero * k(n=n, l=l ,r=r))

def fReal(*, z_list, EA, n, l, r):
    return fReal1(z_list=z_list, EA=EA, n=n, l=l, r=r) + fReal2(z_list=z_list, EA=EA, n=n, l=l, r=r)

def fImag(*, z_list, EA, n, l, r):
    return fImag1(z_list=z_list, EA=EA, n=n, l=l, r=r) - fImag2(z_list=z_list, EA=EA, n=n, l=l, r=r)

def integration_real(*, z_list, EA, DH, n, l, r):

    return cumtrapz(overall(z_list=z_list, DH=DH, n=n, l=l, r=r) * fReal(z_list=z_list, EA=EA, n=n, l=l, r=r), z_list)[-1] #integration real part

def integration_imag(*, z_list, EA, DH, n, l, r):

    return cumtrapz(overall(z_list=z_list, DH=DH, n=n, l=l, r=r) * fImag(z_list=z_list, EA=EA, n=n, l=l, r=r), z_list)[-1] #integration imaginary part

def calculate_array(*, EA_array, DH_array, N, L, R, fR_array, fI_array, F_array):

    for ind_ea, EA in enumerate (EA_array):

        for ind_dh, DH in enumerate (DH_array):              #For each definite value of DH
            z_list= np.linspace(-DH, DH, 100000)    #We define the limit of integration differently,
                                            # the integral is zero outside DH because of the switching function
            result = calculate_single(z_list=z_list, EA=EA, DH=DH, N=N, L=L, R=R)
            F_array[ind_ea][ind_dh] = result

                                                                            #Sum each integral until N
        print(ind_ea)                                      #We set a progress control
        
    return F_array

def calculate_single(*, z_list, EA, DH, N, L, R):

    result = 0

    for n in range(0,N+1):
        for l in range(0, L+1):
            for r in range(0, R+1):

                if n == l == r == 0:
                    continue

                result_integration_real = integration_real(z_list=z_list, EA=EA, DH=DH, n=n, l=l, r=r)

                result_integration_imag = integration_imag(z_list=z_list, EA=EA, DH=DH, n=n, l=l, r=r)
                
                result += cf * (result_integration_real**2 + result_integration_imag**2) #We take the norm times the constant

    
    return result

def function_to_integrate():

    np.cos(np.pi * z_list / (2 * DH)) ** 4 / (n**2 + l**2 + r**2)**0.25 * np.exp(-z_list) / Lzero * (np.cos(z_list * EA) * cos(2 * np.pi * np.exp(-z_list) * (n**2 + l**2 + r**2)**0.5 / Lzero ) + np.sin(z_list * EA) * np.sin(2 * np.pi * exp(-z_list) * (n**2 + l**2 + r**2)**0.5) + 1j * (np.cos(z_list * EA) * np.sin(2 * np.pi * np.exp(-z_list) / Lzero * (n**2 + l**2 + r**2)**0.5) - np.sin(z_list * EA) * np.cos(2 * np.pi * np.exp(-z_list) / Lzero * (n**2 + l**2 + r**2)**0.5) ))

def overall2(z_list, EA, DH, n, l, r):

    return np.cos(np.pi * z_list / (2 * DH)) ** 4 / (n**2 + l**2 + r**2)**0.25 * np.exp(-z_list) / Lzero

def real1(z_list, EA, DH, n, l, r):

    return np.cos(z_list * EA) * np.cos(2 * np.pi * np.exp(-z_list) * (n**2 + l**2 + r**2)**0.5 / Lzero )

def real2(z_list, EA, DH, n, l, r):

    return np.sin(z_list * EA) * np.sin(2 * np.pi * np.exp(-z_list) / Lzero * (n**2 + l**2 + r**2)**0.5)

def imag1(z_list, EA, DH, n, l, r):

    return np.cos(z_list * EA) * np.sin(2 * np.pi * np.exp(-z_list) / Lzero * (n**2 + l**2 + r**2)**0.5)

def imag2(z_list, EA, DH, n, l, r):

    return np.sin(z_list * EA) * np.cos(2 * np.pi * np.exp(-z_list) / Lzero * (n**2 + l**2 + r**2)**0.5)

def integration_real2(z_list, EA, DH, n, l, r):

    return cumtrapz(overall2(z_list, EA, DH, n, l, r) * (real1(z_list, EA, DH, n, l, r) + real2(z_list, EA, DH, n, l, r) ), z_list)[-1]

def integration_imag2(z_list, EA, DH, n, l, r):

    return cumtrapz(overall2(z_list, EA, DH, n, l, r) * (imag1(z_list, EA, DH, n, l, r) - imag2(z_list, EA, DH, n, l, r) ), z_list)[-1]
        
def plot(Y, X, f_array, elev, azim, output_name):

    fig = plt.figure()                              #Define variable of plot
    ax = fig.add_subplot(111, projection='3d')
    
    
    surface = ax.plot_surface(Y, X, f_array , cmap='viridis') #Colour plot
    ax.set_xlabel(r'$E H^{-1}$')                             #Labels
    ax.set_ylabel(r'$\Delta H$')                    #r is for Latex
    ax.set_zlabel('F')
    ax.view_init(elev=elev, azim=azim)                               #Show figure
    plt.savefig(output_name, dpi=650)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='De Sitter')

    # Required arguments
    parser.add_argument('-n', help='N variable', required=True)
    parser.add_argument('-l', help='L variable', required=True)
    parser.add_argument('-r', help='R variable', required=True)

    args = parser.parse_args()

    #sq3 = math.sqrt(3)/(2*H)                  
    N = args.n
    L = args.l
    R = args.r
    
    #Define Grid of parameters
    EA_array = np.linspace(-6,4,50)     #Space of parameters
    DH_array = np.linspace(0.175,3.8,30)    #Is good to take different lenghts so avoid confusion
    X, Y = np.meshgrid(DH_array,EA_array) #Grid definition of parameters
    
    #Canvas
    fR_array = np.zeros( (len(EA_array),len(DH_array)) ) #Define a matrix of zeros = CANVAS (Rows, Columns)
    fI_array = np.zeros( (len(EA_array),len(DH_array)) )
    #ResultadoIntegral = np.zeros( (len(EA_array),len(DH_array)), dtype=np.complex128 ) #We define a matrix of 
                                                            # complex zeros where we gonna join the solution
    F_array = np.zeros( (len(EA_array),len(DH_array)) )       # We gonna use this matrix to turn everything real
    

    F_array = calculate_array(EA_array=EA_array, DH_array=DH_array, N=1, L=1, R=1, fR_array=fR_array, fI_array=fI_array, F_array=F_array)    

    plot(Y=Y, X=X, f_array=F_array, elev=27, azim=+17, output_name="Compact SC deSitter Final Plot 1 DH=(0,3.8,300), EA=(-6, 4, 500), int(100000).png")
    plot(Y=Y, X=X, f_array=F_array, elev=17, azim=-168, output_name="Compact SC deSitter Final Plot 2 DH=(0,3.8,300), EA=(-6, 4, 500), int(100000).png")
    plot(Y=Y, X=X, f_array=F_array, elev=30, azim=+135, output_name="Compact SC deSitter Final Plot 3 DH=(0,3.8,300), EA=(-6, 4, 500), int(100000).png")
    
    plt.show()

    print('Normal termination')

