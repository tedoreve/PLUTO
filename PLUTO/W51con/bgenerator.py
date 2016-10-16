# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:23:47 2016

@author: lingzhi
"""

import numpy as np
import time as ti
from astropy.io import fits
#================================读取密度结构==================================
##双线性插值  
def resize(src,dstsize):#输入src 和size  
    if src.ndim==3:  
        dstsize.append(3)  
    dst=np.array(np.zeros(dstsize),src.dtype)  
    factory=float(np.size(src,0))/dstsize[0]   
    factorx=float(np.size(src,1))/dstsize[1]  
    #print('factory',factory,'factorx',factorx)  
    srcheight=np.size(src,0)  
    srcwidth=np.size(src,1)  
    #print('srcwidth',srcwidth,'srcheight',srcheight)  
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
                t=src[int(fy),int(fx)]*w1+src[int(fy),int(cx)]*w2+src[int(cy),int(fx)]*w3+src[int(cy),int(cx)]*w4
                dst[i,j]=t  
                #print t  
            else:       #映射到原图像刚好是整数的坐标  
                dst[i,j]=(src[int(y),int(x)])  
                # print src[i,j]  
    return dst  
    
def density(filename,width,constant):

    hdulist=fits.open(filename)
    pcdata=hdulist[0].data
    hdulist.close()    
    
    pcdata=resize(pcdata,[width,width])  
    
    pcdata=np.rot90(pcdata,1)
    pcdata=np.flipud(pcdata)
    pcdata=pcdata*constant#13CO转12CO
    return pcdata


#================================加入团块=======================================
def clumps(pcdata,width,number,r,index,c,rho_constant):

    for k in range(number):
        a = np.random.randint(width)
        b = np.random.randint(width)
        for i in range(width):
            for j in range(width):
                c = np.random.randint(r)
                if (i-a)**2+(j-b)**2 < c**2:
                    pcdata[i,j] = c*rho_constant*index**((i-a)**2+(j-b)**2)
    return pcdata

#================================加入星风=======================================
def stellar_wind(pcdata,width,rho_constant,wdir,number,r):
    
    import pyPLUTO as pp
    
    nlinf = pp.nlast_info(w_dir=wdir)    
    D = pp.pload(number,w_dir=wdir)
    print(D.rho.shape)
    
    for i in range(width):
        for j in range(width):
            if (i-256)**2+(j-256)**2<r**2:
                pcdata[i,j] = D.rho[i,j]
    return pcdata
#================================加入磁场=======================================
def magnetism(width,widthi,widthj):

    bx1 = np.zeros([width,width])
    bx2 = np.zeros([width,width])

    for i in range(width):
        for j in range(width):      
            bx1[i,j] = 1e3*(i+widthi)/((i+widthj)**2+(j+widthj)**2)**1.5
            bx2[i,j] = 1e3*(j+widthj)/((i+widthi)**2+(j+widthj)**2)**1.5
    
    bx1  = np.rot90(bx1)
    bx2  = np.rot90(bx2)
    
    bx1 = np.reshape(bx1,width**2,1)
    bx2 = np.reshape(bx2,width**2,1)
    return bx1 , bx2


#================================组合背景=======================================
def combine(components,infilename,outfilename,width,rho_constant,sw,clump,mag):
    
    #密度
    pcdata = np.random.rand(width,width)*rho_constant
    rho    = np.reshape(pcdata,width**2,1) 
    if 'bg' in components:
        pcdata    = density(infilename,width,rho_constant)
        rho       = np.reshape(pcdata,width**2,1)
        total     = rho
    if 'sw' in components:
        wdir,number,r = sw
        pcdata    = stellar_wind(pcdata,width,rho_constant,wdir,number,r)
        rho       = np.reshape(pcdata,width**2,1)
        total     = rho
    if 'clump' in components:
        number,r,index,c = clump
        pcdata    = clumps(pcdata,width,number,r,index,c,rho_constant)
        rho       = np.reshape(pcdata,width**2,1)
        total     = rho
        
    #磁场
    if 'mag' in components:
        widthi,widthj = mag
        bx1 , bx2 = magnetism(width,widthi,widthj)
        total     = np.concatenate((rho,bx1,bx2))
    
    total     = total.astype(float)
    total.tofile(outfilename)
    return pcdata
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
    f.write('1\n')
    f.write('1 0.0 1.0')
    f.close()
    return '空间构造完成！！！'
 
#==============================================================================
if __name__=='__main__':
    print('开始原初构建！！！')
    width = 512
    rho_constant = 20
    sw    = ['../Stellar_Wind/',10,20]    #wdir,number,r
    clump = [2048,4,1.1,0.9]          #number,r,index,c
    mag   = [256,256]                   #widthi,widthj
    combine(['mag','clump'],'W51C.fits','rho0.dbl',width,rho_constant,sw,clump,mag)
    grid('grid0.out',37,512)
    print(ti.asctime())