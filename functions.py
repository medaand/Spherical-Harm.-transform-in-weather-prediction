# -*- coding: utf-8 -*-
"""
Created on Tue Nov  12 09:35:37 2021

@author: Andrea & Marta

In this file we defined the approximations for the Legendre polynomials

"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import math


def c_Abr_Ste(n,m):
    # c coefficients 
    c_n_m = np.sqrt((n-m+1)*(n-m+2)*(n+m+1)*(n+m+2)/((2*n+1)*(2*n+3)**2*(2*n+5)))
    return c_n_m

def d_Abr_Ste(n,m):
    # d coefficients
    d_n_m = (2*n*(n+1)-2*m**2-1)/((2*n-1)*(2*n+3))
    return d_n_m


def P_Abr_Ste(n,m,xk,P_n,P_nminus2):  
# function returning the Legendre polynomial of degree n+2 and order m 
# evaluated at xk, the k-th zero of P_N^m=0
    
# inputs: n = integer, m = integer, xk = real, 
#         P_n = Leg polynomial of degree n and order m evaluated at xk
#         P_nminus2 = Leg polynomial of degree n-2 and order m evaluated at xk
# output: float
    P_n_2 = 1/(c_Abr_Ste(n,m))*((xk**2-d_Abr_Ste(n,m))*P_n-c_Abr_Ste(n-2,m)*P_nminus2)
    return P_n_2
    

def doublefact(n) :
    dfact = 1
    for i in range(1,n+1,2):
        dfact = dfact * i
    return dfact

def dm_P_l(m,l,x):
# gives the m-th derivative of the Legendre polynomial of degree m+l
    if l == 0:
        return doublefact(2*m-1)
    elif l == 1:
        return doublefact(2*m+1)*x
    elif l == 2:
        return doublefact(2*m+1)/2*((2*m+3)*x**2-1)
    elif l == 3:
        return doublefact(2*m+3)/6*x*((2*m+5)*x**2-3)
    
def P_in_vals(m,l,x):
# gives the initial values for the Legendre Polynomials
    P_l_m = np.sqrt(2*l+1)*np.sqrt(math.factorial(l-m)/math.factorial(l+m))*(1-x**2)**(m/2)*dm_P_l(m,l,x)
    return P_l_m


def C_Belusov(n,m):
    # C coefficients 
    C_n_m = np.sqrt(((2*n+1)*(n+m-1)*(m+n-3))/((2*n-3)*(m+n)*(m+n-2)))
    return C_n_m

def D_Belusov(n,m):
    # D coefficients
    D_n_m = np.sqrt(((2*n+1)*(n+m-1)*(n-m+1))/((2*n-1)*(m+n)*(m+n-2)))
    return D_n_m

def E_Belusov(n,m):
    # E coefficients
    E_n_m = np.sqrt(((2*n+1)*(n-m))/((2*n-1)*(m+n)))
    return E_n_m


def P_Belusov(n,m,xk,P_nmin2_mmin2,P_nmin1_mmin2,P_nmin1_m):  
# function returning the Legendre polynomial of degree n+2 and order m 
# evaluated at xk, the k-th zero of P_N^m=0
    
# inputs: n = integer, m = integer, xk = real, 
#         P_n = Leg polynomial of degree n and order m evaluated at xk
#         P_nminus2 = Leg polynomial of degree n-2 and order m evaluated at xk
# output: float
    P_n = C_Belusov(n,m)*P_nmin2_mmin2-D_Belusov(n,m)*xk*P_nmin1_mmin2+E_Belusov(n,m)*xk*P_nmin1_m
    return P_n


def Leg_pol_n(theta,n):
# function giving the ordinary Legendre Poynomial of degree n
# note: equivalent to the associated Legrendre Polynomial with m=0
    M = int(n/2)
    f = 0
    x = np.cos(theta)
    for m in range(0,M+1):
        f = f + (-1)**m*math.factorial(2*n-2*m)/(math.factorial(m)*math.factorial(n-m)*math.factorial(n-2*m))*x**(n-2*m)
    return f/2**n


def dLeg_pol_n(theta,n):
# function giving the derivative of the ordinary Legendre Poynomial of degree n
    M = int(n/2)
    df = 0
    x = np.cos(theta)
    for m in range(0,M+1):
        df = df + (-1)**m*math.factorial(2*n-2*m)/(math.factorial(m)*math.factorial(n-m)*math.factorial(n-2*m))*(n-2*m)*x**(n-2*m-1)
    return df/2**n
 
def d2Leg_pol_n(theta,n):
# function giving the second derivative of the ordinary Legendre Poynomial of degree n
    M = int(n/2)
    df = 0
    x = np.cos(theta)
    for m in range(0,M+1):
        df = df + (-1)**m*math.factorial(2*n-2*m)/(math.factorial(m)*math.factorial(n-m)*math.factorial(n-2*m))*(n-2*m)*(n-2*m-1)*x**(n-2*m-2)
    return df/2**n  
    
    
    
    
    
    
      

    