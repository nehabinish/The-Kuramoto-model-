#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 19:56:56 2020

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt


#FUNCTION TO PLOT AVERAGE OMEGA VS i
def omega_graph(omega_k,N):
    
    '''
    Plots theta values on a circle for time t without N curves 
    
    Parameters:
        ---------------
        omega_k: Velocities of oscillators. NumPy array
        N: No of oscillators
        
        
    See Also:
        ----------------
         np.zeros : Creates a NumPy array with zeroes as values
         
         matplotlib.pyplot.scatter: A scatter plot of y vs. x with varying marker size and/or color.
    
          matplotlib.pyplot.xlabel(xlabel, fontdict=None, labelpad=None, **kwargs)
          Set the label for the x-axis.
            
          matplotlib.pyplot.ylabel(ylabel, fontdict=None, labelpad=None, **kwargs)
          Set the label for the y-axis.
            

     '''
    
    plt.figure()
    omega_av = np.zeros(N)
    for i in range (0,N-1):
         omega_av[i] = np.sum(omega_k[i,:])
         omega_av[i] = 1/N * omega_av[i]
         plt.scatter(i,omega_av[i])
   
    
    plt.xlabel('i')
    plt.ylabel('omega k')
    
    
    
    
# FUNCTION TO PLOT THE ILLUSTRATION OF THE TRANSITION OF CHIMERA STATES
def chimera(z,t):
    
    '''
    Plot illustrating chimera values
   
    Parameters
    ----------
    z: Complex values of order parameter
    t: Time interval for iteration

    
    Returns
    -------
    3-D Graph, 1D Graph
    
    See Also 
    -------
    numpy.meshgrid(*xi, copy=True, sparse=False, indexing='xy')
    Return coordinate matrices from coordinate vectors.
    
    np.arange : Return evenly spaced values within a given interval.
    
    plt.axes():	To create object of the axis
    
    ax.contour3D: To form the contour
    ax.set_xlabel: To label the X-axis
    ax.set_title: To give a title to the plot
    
    view_init(elev=None, azim=None)
    Set the elevation and azimuth of the axes. Rotate the figure
    
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
    circle_angle = np.linspace(0, 2*np.pi, 100)
    radius = np.sqrt(1.0)
    fig, ax = plt.subplots(1)

    x1 = radius*np.cos(circle_angle )
    x2 = radius*np.sin(circle_angle )
    
    plt.plot(x1, x2)
    ax.set_aspect(1)
    plt.grid(linestyle = '--',axis='both')
    plt.xlim(-1.25,1.25)
    plt.ylim(-1.25,1.25)
    
    plt.scatter(z.real,z.imag)
    plt.xlabel('order parameter real')
    plt.ylabel('order parameter imaginary')
    


    plt.figure()
    T,Z = np.meshgrid(t,z)
    ax = plt.axes(projection='3d')
    ax.contour3D(Z.real,Z.imag,T, 50, cmap='viridis')
    ax.set_xlabel('real')
    ax.set_ylabel('imag')
    ax.set_zlabel('time');
    ax.view_init(50,60)
    
  