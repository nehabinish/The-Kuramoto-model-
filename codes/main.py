#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 02:23:18 2020

@author: nehabinish
"""

 #%%
 
import numpy as np

# INITIAL PARAMETERS 
""" 
Parameters
        ----------
        t0 : Start of time interval. (Number)
        tf : End of time interval. (Number)
        step : Total number of values in the interval. (Number)
        N: Total number of oscillators.
        t: Time interval for simulation
        theta: Initialised as NumPy array of zeros
        theta[:,0]: Providing random initial values (between 0-2pi) for first column of theta (time=t0)
        noise: Guassian white noise values between 0 and 1
        w: Providing random natural frequencies between 0 and 1
        M: Number of neighbouring oscillators
        kappa: Critical coupling value
        See Also
        --------
        np.arange : Return evenly spaced values within a given interval.
        np.linspace: Return step numbers of values within a given interva;
        np.random.uniform: Return random float values within a given interval
        np.random.normal: Returns random samples from a normal (Gaussian) distribution
         

"""

t0 = 0
tf = 10
step = 100

N = 20

t = np.linspace(t0,tf,step)

theta=np.zeros((N,len(t)))
theta[:,0] = np.random.uniform(0,2*np.pi,N)

noise = np.random.normal(0,1,size=(N,len(t)))

w = np.random.uniform(-1,1,N)

M = 6



 #%%

# Applying the integrators(Euler,RK2,RK4) to simulate the dynamics of Kuramoto models.
import kuramoto_alpha

kappa = 45

theta_RK4,omega_k = kuramoto_alpha.RK4(kuramoto_alpha.kuramoto,t,theta,N,M,w,kappa)
theta_RK4_mod = theta_RK4%(2*np.pi)
omega_k = omega_k%(2*np.pi)
    
theta_RK2 = kuramoto_alpha.RK2(kuramoto_alpha.kuramoto,t,theta,N,M,w,kappa)
theta_RK2_mod = theta_RK2%(2*np.pi)
        
theta_Euler = kuramoto_alpha.Euler(kuramoto_alpha.kuramoto,t,theta,N,M,w,kappa)
theta_Euler_mod = theta_Euler%(2*np.pi)



 #%%
 
# Plotting the phases against time 

import plot_integrators 
 
plot_integrators.plot_Euler(theta_Euler,N,t)        
plot_integrators.plot_RK2(theta_RK2,N,t)
plot_integrators.plot_RK4(theta_RK4,N,t)

 #%%
 
# Computing and plotting the order parameters

import order_parameter
 
r,phi,z = order_parameter.order(theta_RK4_mod,N,t)
order_parameter.order_graph(theta_RK4_mod,t,r,phi)


 #%%

# Computing and plotting the local Shannon entropy

import shannon_entropy


#Input the index value of the array of time for which the shannon entropy is to be plotted
Time_ = (5,10,20,30,40,50,60,80,99)  
i = 1
      
S = shannon_entropy.Shannon_entropy(theta_RK4_mod,N,t,M)
shannon_entropy.oneD_graph(S,N,t,Time_,i)
shannon_entropy.density_graph(S,N,t,Time_,i)



 #%%
 
 # Plotting the phases on a circle

import plot_phases_circle


#Input time for which the shannon entropy is to be plotted

Time = (t0,tf)

plot_phases_circle.circle_graph(theta_RK4_mod,N,t,Time)

plot_phases_circle.circle_graph_2(theta_RK4_mod,N,t,Time)

plot_phases_circle.density_graph(theta_RK4_mod,N,t)


 #%%

# Illustrating the Chimera states 

import chimera

chimera.omega_graph(omega_k,N)
        
chimera.chimera(z,t)


 #%%
 
# General Case of the Kuramoto model

import kuramoto_modified
'''
Repeat cells by changing input theta to theta_RK4_new after running this cell 
to obtain the values for modified Kuramoto model and to plot graphs.

'''

kappa = 2

K,tau,alpha = kuramoto_modified.random_values(N,t)

theta_RK4_new,omega_k_new = kuramoto_modified.RK4(kuramoto_modified.kuramoto,theta,N,t,w,kappa,M,noise,K,tau,alpha)
theta_RK4_new = theta_RK4_new%(2*np.pi)
omega_k_new = omega_k_new%(2*np.pi)
    
theta_RK2_new = kuramoto_modified.RK2(kuramoto_modified.kuramoto,theta,N,t,w,kappa,M,noise,K,tau,alpha)
theta_RK2_new = theta_RK2_new%(2*np.pi)
        
theta_Euler_new = kuramoto_modified.Euler(kuramoto_modified.kuramoto,theta,N,t,w,kappa,M,noise,K,tau,alpha)
theta_Euler_new = theta_Euler_new%(2*np.pi)




 #%%
 
# Plotting the phases of the Kuramoto model general case against time 

import plot_integrators 
 
plot_integrators.plot_Euler(theta_Euler_new,N,t)        
plot_integrators.plot_RK2(theta_RK2_new,N,t)
plot_integrators.plot_RK4(theta_RK4_new,N,t)


 #%%
