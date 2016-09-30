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
 
##双线性插值  
def resize(src,dstsize):#输入src 和size  
    if src.ndim==3:  
        dstsize.append(3)  
    dst=np.array(np.zeros(dstsize),src.dtype)  
    factory=float(np.size(src,0))/dstsize[0]   
    factorx=float(np.size(src,1))/dstsize[1]  
    print('factory',factory,'factorx',factorx)  
    srcheight=np.size(src,0)  
    srcwidth=np.size(src,1)  
    print('srcwidth',srcwidth,'srcheight',srcheight)  
    for i in range(dstsize[0]):  
        for j in range(dstsize[1]):  
            y=float(i)*factory  
            x=float(j)*factorx  
            if y+1>srcheight:  #越界判断  
                y-=1  
            if x+1>srcwidth:  
                x-=1   
            cy=np.ceil(y)  
            fy=cy-1  
            cx=np.ceil(x)  
            fx=cx-1  
            w1=(cx-x)*(cy-y)  
            w2=(x-fx)*(cy-y)  
            w3=(cx-x)*(y-fy)  
            w4=(x-fx)*(y-fy)      
            if (x-np.floor(x)>1e-6 or y-np.floor(y)>1e-6):   
                t=src[fy,fx]*w1+src[fy,cx]*w2+src[cy,fx]*w3+src[cy,cx]*w4
                dst[i,j]=t  
                #print t  
            else:       #映射到原图像刚好是整数的坐标  
                dst[i,j]=(src[y,x])  
                # print src[i,j]  
    return dst  

print(np.size(pcdata,0))  
pcdata=resize(pcdata,[512,512])  

pcdata=np.rot90(pcdata,1)
pcdata=np.flipud(pcdata)
#plt.imshow(pcdata)
#plt.show()
#==============================================================================

#b=np.fromfile("rho0.dbl",dtype=np.float)
#np.savetxt('xx.txt',b)
width=512
ra=37
rho=np.reshape(pcdata,width**2,1)  #13CO转12CO

x=np.zeros([width,width])
for i in range(width):
    for j in range(width):
        for k in range(width+width):
            if i+j == k:
                x[i,j] = k**1.5/1000000
            
#plt.imshow(x2)
#plt.show()

x1=np.rot90(x,1)+0.008
x2=np.rot90(x,1)+0.008
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
