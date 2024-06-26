# Import packages
import numpy as np
from numpy.fft import rfft
#import matplolib.pyplot as plt

# Question 1: Define function to calculate the fourier coefficients
# Defining the function
def RealToFourier(a, N):
    c = (2/N)*rfft(a)
    return c

def func(x):
    y = np.sin(2*x) - (1/6)*np.cos(4*x)
    return y

# Defining the grid
N = 100
x_i = np.arange(0, 2*np.pi, 2*np.pi/N)

test_func = func(x_i)

test_fourier_coe = RealToFourier(test_func, N)
print('Fourier coefficints:\n', np.around(test_fourier_coe, 3))


# Question 2: Define inverse fourier transform

def FourierToReal(c):