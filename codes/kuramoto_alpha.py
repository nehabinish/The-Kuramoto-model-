#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 18:23:13 2020

@author: nehabinish
"""

import numpy as np
import numba
from numba import jit


# FUNCTION FOR NUMERICAL INTEGRATION USING EULER METHOD

''' 
Numerical integation of the Kuramoto Model

'''

'''
Euler method is a first-order numerical procedure for 
solving ordinary differential equations (ODEs) with a given initial value.

Parameters 
    -------------
    function: kuramoto
    N: Total number of oscillators.
    t: Time interval for simulation
    theta: Initial values of theta
    w: Providing random natural frequencies between 0 and 1
    M: Number of neighbouring oscillators
    kappa: Critical coupling value

'''
@jit(nopython=True)# Imported from Numba module to decrease runtime
def Euler(f,t,theta,N,M,w,kappa):
       
    for n in range(0,len(t)-1):
        delta_t = t[n+1]-t[n]
        theta[:,n+1]=theta[:,n]+f(theta[:,n],t[n],N,t,M,w,kappa)* delta_t
            
    return theta

'''
Returns 
    -------------
    theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
           EULER integration / Phase differences at every point of time.
'''


# FUNCTION FOR NUMERICAL INTEGRATION USING RK2
'''
Error per step for Euler method is very high. 

Runge-Kutta methods are a class of methods which judiciously uses 
the information on the 'slope' at more than one point to extrapolate 
the solution to the future time step

Parameters 
    -------------
    function: kuramoto
    N: Total number of oscillators.
    t: Time interval for simulation
    theta: Initial values of theta
    w: Providing random natural frequencies between 0 and 1
    M: Number of neighbouring oscillators
    kappa: Critical coupling value
'''
@jit(nopython=True) # Imported from Numba module to decrease runtime
def RK2(f,t,theta,N,M,w,kappa):
    
    for n in range(0,len(t)-1):
        delta_t=t[n+1]-t[n]
        theta[:,n+1]=theta[:,n]+f(theta[:,n]+(f(theta[:,n],t[n],N,t,M,w,kappa)*delta_t/2),(t[n]+delta_t/2),N,t,M,w,kappa)*(delta_t)
    
    return theta
'''
Returns 
    -------------
    theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
           RK2 integration / Phase differences at every point of time.
'''



#FUNCTION FOR NUMERICAL INTEGRATION USING RK4

'''
The fourth order Runge-Kutta (RK4) method is more accurate than the RK2 

Parameters 
    -------------
    function: kuramoto
    N: Total number of oscillators.
    t: Time interval for simulation
    theta: Initial values of theta
    w: Providing random natural frequencies between 0 and 1
    M: Number of neighbouring oscillators
    kappa: Critical coupling value
    
See also:
-----------------
   np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)
   
'''

@jit(nopython=True) # Imported from Numba module to decrease runtime
def RK4(f,t,theta,N,M,w,kappa):

    omega_k = np.zeros((N,len(t)))
    
    for n in range(0,len(t)-1):
        delta_t= t[n+1]-t[n]
        k1 = f(theta[:,n],t[n],N,t,M,w,kappa)
        omega_k[:,n] = k1
        k2 = f(theta[:,n]+ k1*delta_t/2,t[n]+ delta_t/2,N,t,M,w,kappa)
        k3 = f(theta[:,n]+ k2*delta_t/2,t[n]+ delta_t/2,N,t,M,w,kappa)
        k4 = f(theta[:,n]+ k3*delta_t,t[n]+ delta_t,N,t,M,w,kappa)
        theta[:,n+1]=theta[:,n]+(k1+2*k2+2*k3+k4) * (delta_t/6)
            

    return theta, omega_k 

'''
Returns 
    -------------
    
    theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
           RK4 integration / Phase differences at every point of time.
    omega_k:  First order derivatives of the phase differences or 
                the velocity of the ith oscillators at every given time t
'''


# FUNCTION TO NUMERICALLY REPRESENT THE KURAMOTO MODEL
'''
Function that numerically represents the Kuramoto Model configuration on which 
the numerical integration is performed.

    Parameters:
        -----------------
        
        theta: Theta value at the time of the working of the integrator
        time: Time value for the working of the integrator
        N: Total number of oscillators.
        t: Time interval for simulation
        theta: Initial values of theta
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        

    See also:
        -----------------
       np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)

'''
@jit(nopython=True) # Imported from Numba module to decrease runtime
def kuramoto(theta,time,N,t,M,w,kappa):
        
        K = np.zeros((N,N))
        sum_ = 0
        theta_dot = np.zeros(N)
        
        for i in range(0,N):
            for j in range (0,N):
                for p in range (1,M): # checking for neighbouring oscillators
                    if(j==i+p or j==i-p): #considering only the nearest neighbour couplings, open chain.
                        K[i][j]=kappa
                        sum_ +=  K[i][j]*np.sin(theta[j]-theta[i])                     
                    else:
                        K[i][j]=0
                        sum_ +=  K[i][j]*np.sin(theta[j]-theta[i])  
                   
            theta_dot[i] = w[i] + (1/N)*sum_
   
        return theta_dot
'''
Returns 
    -------------
    theta_dot : NumPy array with the first order derivative of the Input theta values
    
'''

      

