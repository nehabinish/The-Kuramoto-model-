#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 20:44:44 2020

@author: nehabinish
"""


import numpy as np
import matplotlib.pyplot as plt


# FUNCTION TO COMPUTE ORDER PARAMETERS 'R' AND 'PHI'

def order(theta,N,t):

    """ 
   Computes the Magninute R and angle phi of the complex average of phase oscillators
   
    Parameters
    ----------
    theta : ndarray of shape (N,len(t))
    N: Total number of oscillators.
    t: Time interval for simulation
    
    Returns
    -------
    r: NumPy array of magnitude R
    phi: NumPy array of average phase 
    
    See Also 
    -------
     np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)
     
     numpy.angle(z, deg=False)
     Return the angle of the complex argument.
     
     numpy.absolute: Calculate the absolute value element-wise.
  
    """

    r = np.zeros(len(t))
    phi = np.zeros(len(t))
    z= np.zeros(len(t),dtype = complex)
    sum_order = 0
        
    for n in range(0,len(t)):
        sum_order = np.sum(np.exp(theta[:,n]*1j))
        eq = (1/N)*sum_order
        z[n]=eq
        r[n] = np.absolute(eq)#R(t) given by the absolute value 
        phi[n] = np.angle(eq) #phi(t) given angle of the complex sum value.
     
    print(type(r),type(phi),type(z)) #Checking the type of the return parameters.
    return r,phi,z    





# FUNCION TO PLOT ORDER PARAMETERS t-->r(t) and t-->\phi(t)

def order_graph(theta,t,r,phi):
    
    '''
    Function to plot the order parameters against time 
    
     Parameters
        ----------
        theta : ndarray of shape (N,len(t))
        t: Time interval for simulation
        r: NumPy array of magnitude R
        phi: NumPy array of average phase 
        
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
    
    plt.figure()
    plt.plot(t,r)
    plt.grid(linestyle = '--',axis='both')
    plt.ylim(0,1.2)
    plt.xlabel('time')
    plt.ylabel('R(t)') 

    
    plt.figure()
    plt.plot(t,phi)
    plt.grid(linestyle = '--',axis='both')
    plt.xlabel('time')
    plt.ylabel('\u03C6(t)') 
    
    plt.figure()
    
    plt.plot(t,r,label = 'R(t)')
    plt.grid(linestyle = '--',axis='both')
    plt.plot(t,phi,label='\u03C6(t)')
    plt.grid(linestyle = '--',axis='both')
    plt.xlabel('time')
    plt.ylabel('order parameters') 
