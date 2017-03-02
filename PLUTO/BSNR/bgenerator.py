# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:23:47 2016

@author: lingzhi
"""

import numpy as np



#================================读取密度结构==================================
#import matplotlib.pyplot as plt




#==============================================================================
def magnetism(width):

    bx  = np.zeros([width,width,width])
    x   = np.zeros([width,width,width])


    for i in range(width):
        print(i)
        for j in range(width):
            for k in range(width):
                for l in range(3*width):
                    if i+j+k == l:
                        x[i,j,k] = l

    bx   = np.rot90(x)
    bx   = bx.T/500000
    bx = np.reshape(bx,width**3,1) + 0.001
    
    return bx
#==============================================================================
width=128
ra=37
index=2.4
rho_constant=20

rho    = 1/np.random.power(index,[width,width,width])*rho_constant
rho    = np.reshape(rho,width**3,1)
bx     = magnetism(width)
 
total=np.concatenate((rho,bx,bx,bx))
total=total.astype(float)
total.tofile("rho0.dbl")

#---------------------------------------------------------------------

b=np.linspace(-ra,ra,width+1)
c=np.linspace(1,width,width)
d=np.array([c,b]).T
f=open('grid0.out','w')
f.write('# GEOMETRY:   CARTESIAN\n')
f.write(str(len(b)-1)+'\n')
for i in range(len(c)):
    f.write(str(int(c[i]))+'  '+str(b[i])+'  '+str(b[i+1])+'\n')
f.write(str(len(b)-1)+'\n')
for i in range(len(c)):
    f.write(str(int(c[i]))+'  '+str(b[i])+'  '+str(b[i+1])+'\n')
f.write(str(len(b)-1)+'\n')
for i in range(len(c)):
    f.write(str(int(c[i]))+'  '+str(b[i])+'  '+str(b[i+1])+'\n')
#f.write('1\n')
#f.write('1 1.0 1.0')
f.close()
print(d)
#---------------------------------------------------------------------
