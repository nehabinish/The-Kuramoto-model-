#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:54:50 2021

@author: nehabinish
"""

import numpy as np
import matplotlib.pyplot as plt

plt.figure()
data = np.loadtxt('position.txt')
plt.hist(data, bins=100)
plt.grid()
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()

I = np.average(data)
I = np.average(np.sqrt(data))

print('Integral = ', I)
 
maxl = 200
x = data

N = len(x)       # number of values in array x
m = np.mean(x)
C = np.dot(x-m,x-m)/N

ACF = [1 if l==0 else np.dot(x[:-l] - m, x[l:] - m)/(N-l)/C for l in range(maxl) ]


plt.figure(figsize=(20,10))
plt.plot(ACF)
plt.grid()
plt.xlabel('l')
plt.ylabel('C(l)')


X = np.linspace(0,maxl,maxl)
Y = np.exp(-X/75)
plt.plot(X,Y, ls = ':')
plt.show()

'''
ACF = np.zeros(1000)
l = []

for i in range(0,1000):
    
    l.append(i)
    if i==0 :
        ACF[i] = 1
    else :
        ACF[i] = np.dot(x[:-i]-m,x[i:]-m)/(N-i)/C 
'''        
