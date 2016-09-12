# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:23:47 2016

@author: lingzhi
"""

import numpy as np



#================================读取密度结构==================================
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
hdulist=fits.open('W51C.fits')
pchead=hdulist[0].header
pcdata=hdulist[0].data
hdulist.close()
w=wcs.WCS(pchead)

l1=109.5
l2=108.5
b1=-1.5
b2=-0.5

pcdata=np.rot90(pcdata,1)
pcdata=np.flipud(pcdata)
#plt.imshow(pcdata)
#plt.show()
#==============================================================================

#b=np.fromfile("rho0.dbl",dtype=np.float)
#np.savetxt('xx.txt',b)
width=162
ra=37
rho=np.reshape(pcdata,width**2,1)*45  #13CO转12CO

x=np.zeros([width,width])
for i in range(width):
    for j in range(width):
        for k in range(width+width):
            if i+j == k:
                x[i,j] = k**1.5/10000
            
#plt.imshow(x2)
#plt.show()

x1=np.rot90(x,1)+0.8
x2=np.rot90(x,1)+0.8
bx1=np.reshape(x1,width**2,1)  #13CO转12CO
bx2=np.reshape(x2,width**2,1) 
total=np.concatenate((rho,bx1,bx2))
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
f.write('1\n')
f.write('1 0.0 1.0')
f.close()
print(d)
#---------------------------------------------------------------------
