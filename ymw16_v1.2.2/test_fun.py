import numpy as np                           ## library to work with arrays
import matplotlib.pyplot as plt              ## plotting library
import ctypes
from ctypes import cdll, c_double
lib = cdll.LoadLibrary('./libymw16.so')    ## use our library written in C
#grid = np.linspace (0.01, 1.0, 50)           ## prepare a grid with stepsizes
res_list=[]
lib.d_to_dm.restype = ctypes.c_double          ## Declare a type of return argument (important)

#for i in range (0, len(grid)):               ## for each step size we call our function
#        step = c_double(grid[i])
#        res = lib.value(step)
#        res_list.append(res)

gb   = c_double (0.0)
gl   = c_double (0.0)
dist = c_double (500.0)

print '------------------------'
print lib.d_to_dm(gl, gb, dist)
print '------------------------'
print lib.d_to_dm(gl, gb, dist)
print '------------------------'

print 'Another direction'
gb   = c_double (0.0)
gl   = c_double (90.0)

print lib.d_to_dm(gl, gb, dist)
print lib.d_to_dm(gl, gb, dist)

#plt.plot(grid, res_list)                        ## plot the resylt of numerical integration
#plt.xlabel('Numerical integration step size')
#plt.ylabel('Result of integration')
#plt.savefig('plot.pdf')
#plt.show()
