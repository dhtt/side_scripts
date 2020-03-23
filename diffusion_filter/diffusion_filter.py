#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:51:19 2018

@author: dohoangthutrang
"""

import pandas as pd
import seaborn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


reader = pd.read_csv('noise.csv', delimiter = ',', header = None, index_col = None)
noise = np.array(reader, dtype = float)
#seaborn.heatmap(noise, cmap = cm.magma, cbar = False)

nx = np.shape(noise)[0]
ny = np.shape(noise)[1]

x = np.linspace(0,10, nx)
y = np.linspace(0,10, ny)

X, Y = np.meshgrid(x, y)

def diffusion_function(data, dt, spacing, times):
   D = 1
   k = D*dt/spacing**2
   
   time = 0
   while time < (times+1):
       c = data.copy()
       for i in range(0, nx):
           for j in range(0, ny):
               if i == 0:
                   if j == 0:
                       data[i,j] = k*(c[1,j] + c[0,j] + c[i,1] + c[i,0]) + (1-4*k*c[i,j])
                   elif j == ny-1:
                       data[i,j] = k*(c[1,j] + c[0,j] + c[i,j] + c[i,j-1]) + (1-4*k*c[i,j])
                   else:
                       data[i,j] = k*(c[1,j] + c[0,j] + c[i,j+1] + c[i,j-1]) + (1-4*k*c[i,j])
               elif i == nx - 1:
                   if j == ny - 1:
                       data[i,j] = k*(c[i,j] + c[i-1,j] + c[i,j] + c[i,j-1]) + (1-4*k*c[i,j])
                   elif j == 0:
                       data[i,j] = k*(c[i,j] + c[i-1,j] + c[i,j+1] + c[i,j]) + (1-4*k*c[i,j])
                   else: 
                       data[i,j] = k*(c[i,j] + c[i-1,j] + c[i,j+1] + c[i,j-1]) + (1-4*k*c[i,j])
               
   
               if j == 0:
                   if (0<i<nx-1):
                       data[i,j] = k*(c[i+1,j] + c[i-1,j] + c[i,j+1] + c[i,0]) + (1-4*k*c[i,j])
               elif j == ny - 1:
                   if (0<i<nx-1):
                       data[i,j] = k*(c[i+1,j] + c[i-1,j] + c[i,j] + c[i,j-1]) + (1-4*k*c[i,j])
                       
               if (0 < i < nx-1) and (0 < j < ny-1):
                   noise[i,j] = k*(c[i+1,j] + c[i-1,j] + c[i,j+1] + c[i,j-1]) + (1-4*k*c[i,j])
       time += 1
   
   plt.figure()
   seaborn.heatmap(data, cmap = cm.magma, cbar = False) 
   filename= '2d '+str(dt)+' '+str(times)+'.png'
   plt.savefig(filename, dpi= 96)
   return(data)

diffusion_function(noise, 0.4, 1, 120)
times = [0, 10, 60]
for n in times:
   diffusion_function(noise, 0.10, 1, n)
for n in times:
   diffusion_function(noise, 0.25, 1, n)
for n in times:
   diffusion_function(noise, 0.40, 1, n)
x = 15
y=45
res = x if x<y else y
print(res)