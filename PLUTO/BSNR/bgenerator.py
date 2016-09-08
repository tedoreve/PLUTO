# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:23:47 2016

@author: lingzhi
"""

import numpy as np



#================================读取密度结构==================================
import matplotlib.pyplot as plt




#==============================================================================

#b=np.fromfile("rho0.dbl",dtype=np.float)
#np.savetxt('xx.txt',b)
width=128
ra=30

x=np.zeros([width,width,width])
for i in range(width):
    for j in range(width):
        for k in range(width):
            for l in range(3*width):
                if i+j+k == l:
                    x[i,j,k] = k/10000
            
#plt.imshow(x2)
#plt.show()

x0=np.rot90(x,1)+0.001
bx=np.reshape(x0,width**3,1)  #13CO转12CO
 
total=np.concatenate((bx,bx,bx))
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
