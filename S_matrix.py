# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:16:35 2021

@author: Andrea
"""

import numpy as np
# from scipy.integrate import odeint
# import matplotlib.pyplot as plt
# import pandas as pd
# import math
import functions as fnt
import TP_Newton

N = 50  # horizontal resolution
M = N    # truncating wavenumber 


theta = np.linspace(0,2*np.pi,1000)

in_guess = np.linspace(np.pi+0.1,2*np.pi-0.1,N)   # vector for equidistant initial guesses
itmax = 200    # /!\ to be fixed
tol = 10**-8   # /!\ to be fixed
xk = np.zeros([N,1])



for i in range(0,len(in_guess)):
    # calculate and store in an array of the roots P_N
    xk[i] = TP_Newton.newton(fnt.Leg_pol_n,fnt.dLeg_pol_n,in_guess[i],tol,itmax,N) 
    
# print(xk)

S = np.zeros([N,N])  # initialization of the matrix

S_store = []         # initializations of the tuple containing all matrixes

# building the matrix of order m=0
for i in range(0,N):
    for j in range(1,N+1):
        S[i][j-1] = fnt.Leg_pol_n(xk[i],j)
        
S_store.append(S)  # storing the matrix of order m=0

S = np.zeros([N,N])  # initialization of the matrix

# building the matrix of order m =1
for i in range(0,N):
    for j in range(1,N+1):
        S[i][j-1] = fnt.dLeg_pol_n(xk[i],j)*np.sqrt(1-np.cos(xk[i])**2)
        
S_store.append(S)  # storing the matrix of order m=1

S = np.zeros([N,N])  # initialization of the matrix

# building the matrx of order m=2
for i in range(0,N):
    for j in range(2,N+1):
        if j == 2:
            S[i][j-1] = fnt.d2Leg_pol_n(xk[i],j)*(1-np.cos(xk[i])**2)
            
        else:
            S[i][j-1] = fnt.P_Belusov(j,2,xk[i],S_store[0][i][j-3],S_store[0][i][j-2],S[i][j-2])
            
S_store.append(S) #storing


# building and storing all matrixes with Belusov formulation (7)
for l in range(3,M+1):
    S = np.zeros([N,N])
    for i in range(0,N):
        for j in range(l,N+1):
            S[i][j-1] = fnt.P_Belusov(j,l,xk[i],S_store[l-2][i][j-3],S_store[l-2][i][j-2],S[i][j-2])
    S_store.append(S)
    
    
    
    
    
    
    
    