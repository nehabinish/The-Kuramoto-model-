#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 16:31:56 2020

@author: nehabinish
"""


import numpy as np
import matplotlib.pyplot as plt




# FUNCTION TO PLOT THE N-CURVES t --> theta_i(t)
def N_curves(theta,N,t):
    
    '''
    Plots the N curves on the same graph as the circle.
    
    Parameters:
        ---------------
         N: Total number of oscillators.
         t: Time interval for simulation
         theta: NumPy array of theta obtained from Euler, RK2 and RK4 respectively.
         
    See Also:
        ----------------
         np.arange : Return evenly spaced values within a given interval.
         
         matplotlib.pyplot.plot(*args, scalex=True, scaley=True, data=None, **kwargs)
         Plot y versus x as lines and/or markers.
         
    '''
    
    oscillators = np.arange(0,N,1)
    
    for i in oscillators: 
        
        plt.plot(t,theta[i,:])



# FUNCTION TO REPRESENT GRAPH ON A CIRCLE AT TIME 't'
def circle_graph(theta,N,t,Time):
    '''
    
    Plots theta values on a circle for time t with the N curves
    
    Parameters:
        ---------------
         N: Total number of oscillators.
         t: Time interval for simulation
         theta: NumPy array of theta obtained from Euler, RK2 and RK4 respectively.
         Time: Given time for which theta values are plotted.
         
    See Also:
        ----------------
         np.arange : Return evenly spaced values within a given interval.
         
         numpy.where(condition[, x, y])
         Return elements chosen from x or y depending on condition.
         
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
    for n in Time:
        
        x=np.where(t==n)[0]

        circle_angle = np.linspace(0, 2*np.pi, 100)
        radius = np.sqrt(4*(np.pi)**2)
        plt.figure()
        fig, ax = plt.subplots(1)

        x1 = radius*np.cos(circle_angle )
        x2 = radius*np.sin(circle_angle)

        plt.plot(x1, x2)
        ax.set_aspect(1)
        plt.grid(linestyle = '--',axis='both') #Drawing a circle with radius 2pi
     
        plt.plot(radius*np.cos(theta[:,x]),radius*np.sin(theta[:,x]),'*') #Drawing phases on the circle
        plt.title("time = {}".format(n))
        
        #Calling the function to plot N_curves in the same graph
        N_curves(theta,N,t)
 
    

def circle_graph_2(theta,N,t,Time):
    
    '''
    
    Plots theta values on a circle for time t without N curves 
    
    Parameters:
        ---------------
         N: Total number of oscillators.
         t: Time interval for simulation
         theta: NumPy array of theta obtained from Euler, RK2 and RK4 respectively.
         Time: Given time for which theta values are plotted.
         
    See Also:
        ----------------
         np.arange : Return evenly spaced values within a given interval.
         
         numpy.where(condition[, x, y])
         Return elements chosen from x or y depending on condition.
         
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
    
    oscillators = np.arange(0,N,1)
    
    for n in Time:
        
        x=np.where(t==n)[0]
        
    
        circle_angle = np.linspace(0, 2*np.pi, 100)
        radius = np.sqrt(4*(np.pi)**2)
        plt.figure()
        fig, ax = plt.subplots(1)

        x1 = radius*np.cos(circle_angle )
        x2 = radius*np.sin(circle_angle )

        plt.plot(x1, x2)
        ax.set_aspect(1)
        plt.grid(linestyle = '--',axis='both')
     
        plt.plot(radius*np.cos(theta[:,x]),radius*np.sin(theta[:,x]),'*')  
        plt.title("time = {}".format(n))
        
        for a in oscillators:
            plt.plot([0, radius*np.cos(theta[a][x])], [0, radius*np.sin(theta[a][x])])
            plt.title("time = {}".format(n))



# FUNCTION TO DRAW THE DENSITY GRAPH (i,t) --> theta_i(t)
def density_graph(theta,N,t):
    
    """ 
    Density graph to plot (i,t) ---> theta(i,t) at a given time
   
    Parameters
    ----------
    theta: NumPy array with phase differences
    N:  The  total  number  of  oscillators
    t: Time interval for iteration of the Kuramoto model

    
    Returns
    -------
    3-D Graph
    
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
    
 
    """
    i =np.arange(0,N,1)
    T,I = np.meshgrid(t,i)
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(T,I,theta, 50, cmap='viridis')
    ax.set_xlabel('time')
    ax.set_ylabel('i')
    ax.set_zlabel('$\u03B8^i(t)$');
    ax.view_init(25,60)



