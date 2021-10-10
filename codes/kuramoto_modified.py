#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 20:37:31 2020

@author: nehabinish
"""

import numpy as np
import numba
from numba import jit


'''
Applying the integrators(Euler,RK2,RK4) to simulate the dynamics of Kuramoto 
models.

'''

# FUNCTION FOR NUMERICAL INTEGRATION USING EULER METHOD
@jit(nopython=True)# Imported from Numba module to decrease runtime
def Euler(f,theta,N,t,w,kappa,M,noise,K,tau,alpha):
    
    '''
    Euler method is a first-order numerical procedure for 
    solving ordinary differential equations (ODEs) with a given initial value.

    Parameters 
        -------------
        function kuramoto
        theta: Initial values of theta
        N: Total number of oscillators.
        t: Time interval for simulation   
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        noise: Guassian white noise
        K: Coupling Matrix
        tau: Delay Matrix
        alpha: Dephasing Matrix
        
    
    Returns 
        -------------
        theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
              Euler integration / Phase differences at every point of time.


    '''
       
    for n in range(0,len(t)-1):
        delta_t = t[n+1]-t[n]
        theta[:,n+1]=theta[:,n]+f(theta[:,n],t[n],noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)* delta_t
            
    return theta




# FUNCTION FOR NUMERICAL INTEGRATION USING RK2
@jit(nopython=True)# Imported from Numba module to decrease runtime
def RK2(f,theta,N,t,w,kappa,M,noise,K,tau,alpha):
    
    '''
    Runge-Kutta methods are a class of methods which judiciously uses 
    the information on the 'slope' at more than one point to extrapolate 
    the solution to the future time step

    Parameters 
        -------------
        function kuramoto
        theta: Initial values of theta
        N: Total number of oscillators.
        t: Time interval for simulation   
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        noise: Guassian white noise
        K: Coupling Matrix
        tau: Delay Matrix
        alpha: Dephasing Matrix
        
    
    Returns 
    -------------
        theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
              RK2 integration / Phase differences at every point of time.

    '''
    
    for n in range(0,len(t)-1):
        delta_t=t[n+1]-t[n]
        theta[:,n+1]=theta[:,n]+f(theta[:,n]+(f(theta[:,n],t[n],noise[:,n],\
                     theta,n,N,M,w,kappa,K,tau,alpha)*delta_t/2),(t[n]+delta_t/2),noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)\
                     *(delta_t)
    
    return theta




#FUNCTION FOR NUMERICAL INTEGRATION USING RK4
@jit(nopython=True)# Imported from Numba module to decrease runtime
def RK4(f,theta,N,t,w,kappa,M,noise,K,tau,alpha):
    
    '''
    More precise order of Runge Kutta integration

    Parameters 
        -------------
        function kuramoto
        theta: Initial values of theta
        N: Total number of oscillators.
        t: Time interval for simulation   
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        noise: Guassian white noise
        K: Coupling Matrix
        tau: Delay Matrix
        alpha: Dephasing Matrix
        
    
    Returns 
    -------------
        theta: Numpy Array - Solution for the solved ODE of the Kuramoto model using 
              RK4 integration / Phase differences at every point of time.
        omega_k:  First order derivatives of the phase differences or 
                the velocity of the ith oscillators at every given time t   
                
     See also:
     -----------------
       np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)
       
    '''

    omega_k = np.zeros((N,len(t)-1))

    for n in range(0,len(t)-1):
        delta_t= t[n+1]-t[n]
        k1 = f(theta[:,n],t[n],noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)
        omega_k[:,n] = k1
        k2 = f(theta[:,n]+ k1*delta_t/2,t[n]+ delta_t/2,noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)
        k3 = f(theta[:,n]+ k2*delta_t/2,t[n]+ delta_t/2,noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)
        k4 = f(theta[:,n]+ k3*delta_t,t[n]+ delta_t,noise[:,n],theta,n,N,M,w,kappa,K,tau,alpha)
        theta[:,n+1]=theta[:,n]+(k1+2*k2+2*k3+k4) * (delta_t/6)
            

    return theta, omega_k



#FUNCTIONS TO BUILD K, τ AND α AS RANDOM MATRICES
@jit(nopython=True)# Imported from Numba module to decrease runtime
def random_values(N,t):
    
    K = np.random.uniform(0,2*np.pi,size=(N,N))
    tau = np.random.randint(0,2*np.pi,size=(N,N))
    alpha = np.random.uniform(0,2*np.pi, size=(N,N))
    
    return K,tau,alpha



# FUNCTION TO NUMERICALLY REPRESENT THE KURAMOTO MODEL
'''
In the previous version, the kuramoto model was represented without noises, dephasing nor delay. 
However in this we add the noise array at time t as well as the matrix of theta values that 
will help us introduce the time delay tau[i][j] and a dephasing matrix alpha[i][j]

Parameters:
        -----------------
        theta_i: Theta value at the time of the working of the integrator
        time: Time value for the working of the integrator
        eta: Guassian white noise
        theta_matrix: NumPy array of phase difference (N,len(t))
        n: Array Index of the time for the working of the integrator
        N: Total number of oscillators.
        theta: Initial values of theta
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        K: Coupling Matrix
        tau: Delay Matrix
        alpha: Dephasing Matrix
        

    See also:
        -----------------
       np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)


'''

@jit(nopython=True)# Imported from Numba module to decrease runtime
def kuramoto(theta_i,time,eta,theta_matrix,n,N,M,w,kappa,K,tau,alpha):
    

        sum_ = 0
        theta_dot = np.zeros(N)
        
        for i in range(0,N-1):
            for j in range (0,N-1):
                for p in range (1,M):
                    theta_j = theta_matrix[:,n-tau[i][j]]
                    if(j==i+p or j==i-p):
                        K[i][j]=kappa
                        sum_ +=  K[i][j]*np.sin(theta_j[j]-theta_i[i]+\
                                                alpha[i][j]) + eta[i]                    
                    else:
                        K[i][j]=0
                        sum_ +=  K[i][j]*np.sin(theta_j[j]-theta_i[i]+\
                                                alpha[i][j]) + eta[i]  
                   
            theta_dot[i] = w[i] + (1/N)*sum_
   
        return theta_dot
    
    



