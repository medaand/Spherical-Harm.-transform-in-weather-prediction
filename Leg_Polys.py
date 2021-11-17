# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 20:17:56 2021

@author: Andrea
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import functions as fnt

xvals = np.linspace(0,2*np.pi,1000)

def legpolyns(x,n):
    leg = sp.legendre(n)
    P_n = leg(np.cos(x))
    return P_n
    

N = 4 # number of polynomials
plt.figure(1)
for i in range(0,N):
    func = legpolyns(xvals, i)
    dfunc =  np.polyder(sp.legendre(i), 1)
    plt.plot(xvals,dfunc(np.cos(xvals)), label = 'n = ' + str(i))
   

plt.title(str(N) + '  Legendre polynomials of the cosine')
plt.grid()
plt.legend(loc='upper left')
plt.show()

plt.figure(2)
for i in range(0,N):
    func1 = fnt.dLeg_pol_n(xvals, i)      # the derivative
    plt.plot(xvals,func1, label = 'n = ' + str(i))
   

plt.title(str(N) + '  Legendre polynomials by us')
plt.grid()
plt.legend(loc='upper left')
plt.show()

xvals = np.linspace(-1,1,1000)

def legpolyns(x,n):
    leg = sp.legendre(n)
    P_n = leg(x)
    return P_n
    

N = 5 # number of polynomials
plt.figure(2)
for i in range(0,N):
    func = legpolyns(xvals, i)
    plt.plot(xvals,func, label = 'n = ' + str(i))
    
plt.title(str(N) + '  Legendre polynomials')
plt.grid()
plt.legend(loc='upper left')
plt.show()