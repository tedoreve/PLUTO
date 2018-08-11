# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:23:47 2016

@author: lingzhi
"""

import numpy as np
import time as ti
from astropy.io import fits

#================================加入星风=======================================
def stellar_wind(wdir,number):

    import pyPLUTO as pp

    pp.nlast_info(w_dir=wdir)
    D = pp.pload(number,w_dir=wdir)
    print(D.rho.shape)

    rho = D.rho
    vx1 = D.vx1
    vx2 = D.vx2

    rho = np.transpose(rho)

    vx1 = np.transpose(vx1)
    vx2 = np.transpose(vx2)

                    
    return rho,vx1,vx2
#================================加入磁场=======================================
def toff(f):
    def wrapper(*args):
        start = ti.time()
        f(*args)
        end   = ti.time()
        print(end-start)
    return wrapper

def f(i,j,k):
    return i*2+j*2

#@toff
def magnetism(width):
    x = np.fromfunction(f,(width,width,width))/500000
    x = np.rot90(x,k=1)
    x = np.transpose(x)
    x = np.reshape(x,width**3,1)
    x = x + 0.001
    return x*0, x, x


#================================组合背景=======================================
def combine(components,infilename,outfilename,width,index,rho_constant,sw,clump,mag):

    if 'sw' in components:
        wdir,number = sw
        rho,vx1,vx2    = stellar_wind(wdir,number)
        rho       = np.reshape(rho,(1,width**2))
        vx1       = np.reshape(vx1,(1,width**2))
        vx2       = np.reshape(vx2,(1,width**2))
        total     = np.concatenate((rho,vx1,vx2))
#
    total     = total.astype(float)
    total.tofile(outfilename)
    return total 
#================================网格定义=======================================
def grid(outfilename,ra,width):

    b=np.linspace(-ra,ra,width+1)
    c=np.linspace(1,width,width)
    f=open(outfilename,'w')
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
#    f.write('1\n')
#    f.write('1 0.0 1.0')
    f.close()
    return '空间构造完成！！！'

#==============================================================================
if __name__=='__main__':
    print('开始原初构建！！！')
    width = 512
    index = 1.0         
    u     = 1.3
    rho_constant = 0.21*u
    sw    = ['./',10]    #wdir,number,r
#    sw    = ['../SW2/',10]
    clump = [200,10,1.0,50.0]             #number,r,index,e
    mag   = 3.2                           #widthi,widthj
    total = combine(['sw'],'W51C.fits','rho0.dbl',width,index,rho_constant,sw,clump,mag)
    grid('grid0.out',0.1,width)
    print(ti.asctime())
    
