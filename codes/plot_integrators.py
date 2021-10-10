#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 20:21:56 2020

@author: nehabinish
"""


import matplotlib.pyplot as plt

'''
Euler method is a first-order numerical procedure for 
solving ordinary differential equations (ODEs) with a given initial value.

Parameters 
    -------------
    N: Total number of oscillators.
    t: Time interval for simulation
    theta: NumPy array of theta obtained from Euler, RK2 and RK4 respectively.
    
    See also
    -------------
    matplotlib.pyplot.plot(*args, scalex=True, scaley=True, data=None, **kwargs)
    Plot y versus x as lines and/or markers.
    
    matplotlib.pyplot.grid(b=None, which='major', axis='both', **kwargs)
    Configure the grid lines.
    
    matplotlib.pyplot.xlabel(xlabel, fontdict=None, labelpad=None, **kwargs)
    Set the label for the x-axis.
    
    matplotlib.pyplot.ylabel(ylabel, fontdict=None, labelpad=None, **kwargs)
    Set the label for the y-axis.
    
    matplotlib.pyplot.title(label, fontdict=None, loc='center', pad=None, **kwargs)
    Set a title for the axes.



'''

def plot_Euler(theta,N,t):
    plt.figure()
    for i in range (0,N-1):
        plt.plot(t,theta[i,:])
        plt.grid(linestyle = '--',axis='both')
        plt.xlabel('time')
        plt.ylabel('\u03B8') #unicode character for theta
        plt.title("Euler")
                

def plot_RK2(theta,N,t):               
    plt.figure()
    for i in range (0,N-1): 
        plt.plot(t,theta[i,:])
        plt.grid(linestyle = '--',axis='both')
        plt.xlabel('time')
        plt.ylabel('\u03B8') #unicode character for theta
        plt.title("RK2")
              

def plot_RK4(theta,N,t):              
    plt.figure()
    for i in range (0,N-1):
        plt.plot(t,theta[i,:])
        plt.grid(linestyle = '--',axis='both')
        plt.xlabel('time')
        plt.ylabel('\u03B8') #unicode character for theta
        plt.title('RK4') 
	plt.show()	




