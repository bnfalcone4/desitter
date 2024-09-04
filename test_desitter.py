#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:35:02 2024

@author: alvaro
"""

from desitter import *

diff = 0.000001


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


def test_calculate_array1():
    
    EA_array = np.linspace(-10,0,50)
    DH_array = np.array([3.5])
    
    #Canvas
    fR_array = np.zeros( (len(EA_array),len(DH_array)) ) #Define a matrix of zeros = CANVAS (Rows, Columns)
    fI_array = np.zeros( (len(EA_array),len(DH_array)) )
    #ResultadoIntegral = np.zeros( (len(EA_array),len(DH_array)), dtype=np.complex128 ) #We define a matrix of 
                                                            # complex zeros where we gonna join the solution
    F_array = np.zeros( (len(EA_array),len(DH_array)) )  
    
    F_array = calculate_array(EA_array=EA_array, DH_array=DH_array, N=1, L=1, R=1, fR_array=fR_array, fI_array=fI_array, F_array=F_array)

    print('F_array', F_array)

    #X, Y = np.meshgrid(DH_array,EA_array)

    #plot(Y, X, F_array, elev=30, azim=20, output_name='test.pdf')


def test_fimag():

    DH = 1

    z_list= np.linspace(-DH, DH, 10)

    fimag = fImag(z_list=z_list, EA=1, n=1, l=1, r=1)

    print(fimag)


def test_fimag2():

    DH = 1

    z_list= np.linspace(-DH, DH, 10)

    fimag2 = fImag2(z_list=z_list, EA=1, n=1, l=1, r=1)

    print(fimag2)


def test_switch():

    DH = 1

    z_list= np.array([0, DH, 2 * DH, 3 * DH, 4 * DH])

    output = switch(z_list=z_list, DH=DH)

    check = np.array([1.0, 0.0, 1.0, 0.0, 1.0])

    for i, j, in zip(output, check):

        assert (i - j) ** 2 < diff


def test_overal():

    DH = 1

    z_list= np.array([DH, 3 * DH])

    output = overall(z_list=z_list, DH=DH, n=1, l=1, r=1)

    check = np.array([0.0, 0.0])

    for i, j, in zip(output, check):

        assert (i - j) ** 2 < diff


def test_array():
    
    EA_array = np.linspace(0, 9, 10)
    DH_array = np.linspace(0, 9, 10)
    
    #Canvas
    fR_array = np.zeros( (len(EA_array),len(DH_array)) ) #Define a matrix of zeros = CANVAS (Rows, Columns)
    fI_array = np.zeros( (len(EA_array),len(DH_array)) )
    #ResultadoIntegral = np.zeros( (len(EA_array),len(DH_array)), dtype=np.complex128 ) #We define a matrix of 
                                                            # complex zeros where we gonna join the solution
    F_array = np.zeros( (len(EA_array),len(DH_array)) )     
    
    F_array[0] = 1

    print(F_array)


def test_overall2():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = overall(z_list=z_list, DH=DH, n=n, l=l, r=r)

    output2 = overall2(z_list, EA, DH, n, l, r)

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_real1():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fReal1(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = real1(z_list, EA, DH, n, l, r)

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_real2():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fReal2(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = real2(z_list, EA, DH, n, l, r)

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_imag1():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fImag1(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = imag1(z_list, EA, DH, n, l, r)    

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_imag2():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fImag2(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = imag2(z_list, EA, DH, n, l, r)

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_real():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fReal(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = real1(z_list, EA, DH, n, l, r) + real2(z_list, EA, DH, n, l, r)

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff   


def test_imag():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = fImag(z_list=z_list, EA=EA, n=n, l=l, r=r)

    output2 = imag1(z_list, EA, DH, n, l, r) - imag2(z_list, EA, DH, n, l, r)  

    for i, j in zip(output1, output2):

        assert (i-j)**2 < diff


def test_loop_real():

    N = 3
    L = 3
    R = 3

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1

    for n in range(N):
        for l in range(L):
            for r in range(R):
                sum_real1 = fReal(z_list=z_list, EA=EA, n=n, l=l, r=r)
                sum_real2 = real1(z_list, EA, DH, n, l, r) + real2(z_list, EA, DH, n, l, r)

                assert np.sum(sum_real1 ** 2 -sum_real2 ** 2) < diff


def test_loop_imag():

    N = 3
    L = 3
    R = 3

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1

    for n in range(N):
        for l in range(L):
            for r in range(R):
                sum_real1 = fImag(z_list=z_list, EA=EA, n=n, l=l, r=r)
                sum_real2 = imag1(z_list, EA, DH, n, l, r) + imag2(z_list, EA, DH, n, l, r)

                assert np.sum(sum_real1 ** 2 -sum_real2 ** 2) < diff


def test_integration_real():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = integration_real(z_list=z_list, EA=EA, DH=DH, n=n, l=l, r=r)
    output2 = integration_real2(z_list, EA, DH, n, l, r)

    print(output1, output2)

    assert (output1 - output2) ** 2 < diff


def test_integration_imag():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    n = 1
    l = 1
    r = 1

    output1 = integration_imag(z_list=z_list, EA=EA, DH=DH, n=n, l=l, r=r)
    output2 = integration_imag2(z_list, EA, DH, n, l, r)

    print(output1, output2)

    assert (output1 - output2) ** 2 < diff


def test_calculate_single():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    DH = 1
    N = 1
    L = 1
    R = 1    

    result = calculate_single(z_list=z_list, EA=EA, DH=DH, N=N, L=L, R=R)

    print(result)


def test_calculate_single2():

    z_list = np.linspace(-10, 10, 1000)

    DH = 1
    N = 1
    L = 1
    R = 1

    for EA in range(-10, 10, 1):

        result = calculate_single(z_list=z_list, EA=EA, DH=DH, N=N, L=L, R=R)

        print(result)

        assert result != 0


def test_calculate_single3():

    z_list = np.linspace(-10, 10, 1000)

    EA = 1
    N = 1
    L = 1
    R = 1

    for DH in range(1, 10, 1):

        result = calculate_single(z_list=z_list, EA=EA, DH=DH, N=N, L=L, R=R)

        print(result)

        assert result != 0.0


def test_diff():

    a = np.linspace(-10, 0, 1000)
    b = np.linspace(0, 10, 1000)

    print(a, b)

    for i, j in zip(a, b):

        assert i != j

    c = np.linspace(-10, 0, 1000)

    for i, j in zip(a, c):

        assert (i-j)**2 < diff


def test_enumerate_array():

    EA_array = np.linspace(-10,0,10)

    print(EA_array)

    for ind_ea, EA in enumerate(EA_array):
        print(ind_ea, EA)



