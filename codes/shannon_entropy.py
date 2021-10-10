#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 18:37:42 2020

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt
import numba
from numba import jit


# Function to calculate the entropy.
jit(nopython=True) # To enhance the speed and decrease runtime 
def Shannon_entropy(theta_i,N,t,M):
    
    """ 
    Computes the local Shannon entropy for all the oscillators
   
    Parameters
    ----------
    theta : Phase difference. ndarray of shape (N,len(t))
    N: Total number of oscillators.
    t: Time interval for simulation
    M: Total number of neighbouring oscillators.
    
    Returns
    -------
    S: NumPy array containing the local Shannon entropy of all the oscillators 
    
    See Also 
    -------
     np.zeros: Return a new array of given shape and type, filled with zeros(used for intilisation)
     
     numpy.floor: Return the floor of the input, element-wise.
                    The floor of the scalar x is the largest integer i, 
                    such that i <= x. It is often denoted as \lfloor x \rfloor.
                    
     
     astype(data_type): to change the data type of a numpy array
     
     numpy.absolute: Calculate the absolute value element-wise.
  
    """
    
    q=np.zeros((N,len(t)))
    q[:,0] = (2*M)+1 #Setting intial maximum elongation as the total reach of an oscillator
    
    
    for n in range (1,len(t)-1): #Finding the rest of the q values
        q[:,n] = q[:,0]*np.cos(theta_i[:,n])


    q = np.floor(q) #Rounding the elongation value
    q = q.astype(int) #Converting float values to int so it can be provided in a range
    q = np.absolute(q) # Finding absolute values of q
    q[q==0]=2 #Setting the minimun elongation of q as 2
    
    print("q")
    print(q)
    
    S = np.zeros((N,len(t)))
    
    a = 0
    

    theta_k = np.zeros(N)
    
    for n in range(0,len(t)):
        for i in range (0,N):
            
            for k in np.arange(i-M,i+M,1): #Checking the neighbouring oscillators, Here, closed chain
                
                if (i+M > N): # Modulo N since it's a closed chain.
                    k = (i+M)%N
                           
                    if (theta_i[k][n] >=(2*np.pi*a/q[i][n])) and (theta_i[k][n]<= (2*np.pi*(a+1)/q[i][n])):
                        theta_k[k] = theta_i[k][n]
                    
                else:
                    #Checking if theta_k falls in this set of values
                    if ((2*np.pi*a/q[i][n]) <= (theta_i[k][n])) and ((theta_i[k][n])<= (2*np.pi*(a+1)/q[i][n])):
                        theta_k[k]=theta_i[k][n]
                        
            p_a = len(theta_k)/((2*M)+1) #Calculating probability to be the fraction of oscillators with phase in range(2*np.pi*a/q[i][n],2*np.pi*(a+1)/q[i][n])

            
            for p in range(a,q[i][n]-1):
                S[i][n] += - p_a * np.log(p_a)  
                                  
       
    return S
   
    

# Function to draw 1-D graph for Shannon entropy

def oneD_graph(S,N,t,Time_,i):
    
        
    """ 
    1-D graph to plot Shannon entropy at a given time
   
    Parameters
    ----------
    S: The Local Shannon entropy values.Numpy array
    N:  The  total  number  of  oscillators
    t: Time interval for iteration of the Kuramoto model
    Time:  A set of time for which the Shannon entropy is to be plotted.
    i: ith oscillator for which the entropy change is to be plotted
    
    Returns
    -------
    1-D Graph
    
    See Also 
    -------
    matplotlib.pyplot.scatter: A scatter plot of y vs. x with varying marker size and/or color.
 
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
    
    np.arange : Return evenly spaced values within a given interval.
    
    numpy.where(condition[, x, y])
         Return elements chosen from x or y depending on condition.
  
    """
    
    oscillators = np.arange(0,N,1)
  
    plt.figure()
    
    for x in Time_:
         
        plt.figure()
        
        for i in oscillators:
            plt.scatter(i,S[i][x])
            plt.title("t_index={}".format(x))
            plt.xlabel("i")
            plt.ylabel("$S_{i}^{q,n}(t)$")   
    
       
    plt.figure()

    plt.plot(t,S[i,:])
    plt.title("Shannon entropy")
    plt.xlabel("t")
    plt.ylabel("$S_{i}^{q,n}(t)$") 


# Function to draw 3-D graph for Shannon entropy
def density_graph(S,N,t,Time_,i):
    
    """ 
    Density graph to plot (i,t) ---> Shannon entropy(i,t) at a given time
   
    Parameters
    ----------
    S: The Local Shannon entropy values.Numpy array
    N:  The  total  number  of  oscillators
    t: Time interval for iteration of the Kuramoto model
    Time:  A set of time for which the Shannon entropy is to be plotted.
    i: ith oscillator for which the entropy change is to be plotted
    
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
    ax.contour3D(T,I,S, 50, cmap='viridis')
    ax.set_xlabel('time')
    ax.set_ylabel('i')
    ax.set_zlabel('$S^i(t)$');
    ax.view_init(25,60)
