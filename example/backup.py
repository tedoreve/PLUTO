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
hdulist=fits.open('CTB109.fits')
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
plt.imshow(pcdata)#,origin='lower',interpolation='nearest',extent=[l1,l2,b1,b2])
#==============================================================================

#b=np.fromfile("rho0.dbl",dtype=np.float)
#np.savetxt('xx.txt',b) 
a=np.reshape(pcdata,40000,1)
a=a.astype(float)
aa=np.random.uniform(0,1,size=40000)
a.tofile("rho0.dbl")

#---------------------------------------------------------------------

b=np.linspace(-0.5,0.5,201)
c=np.linspace(1,200,len(b)-1)
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

#---------------------------------------------------------------------
#Spectral_Index=-5/3
#Minimum_Lambda = 4
#Norm      = 0.1
#IR_Cutoff = 2.0
#NXLim     = 2*(int(len(c)/2/Minimum_Lambda))
#NYLim     = 2*(int(len(c)/2/Minimum_Lambda))
#NXLim_2   = NXLim/2
#NYLim_2   = NYLim/2
#
#Amplitude = np.zeros([NXLim + 1, NYLim + 1])
#Phase     = np.zeros([NXLim + 1, NYLim + 1])
#twoPI     = np.pi*2
#x_end     = 1
#y_end     = 1
#for i in range(NXLim):
#    for j in range(NYLim):
#        x1r  = np.random.uniform(0,1,size=1)
#        x2r  = np.random.uniform(0,1,size=1)
#        xdbl = np.sqrt(-2*np.log(x1r))*np.cos(twoPI*x2r)
#        Phase[i,j]=twoPI*np.random.uniform(0,1,size=1)
#        scrh1= IR_Cutoff**2+(twoPI*(i-NXLim_2)/x_end)**2+(twoPI*(i-NYLim_2)/y_end)**2
#        Amplitude[i,j]=np.sqrt(Norm*scrh**(Spectral_Index/2))*xdbl
#for j2 in range():












