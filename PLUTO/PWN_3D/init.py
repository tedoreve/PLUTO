# -*- coding: utf-8 -*-
"""
Created on Tue May  9 20:09:42 2017

@author: tedoreve
"""

from astropy import constants as con
from astropy import units as un
import math as ma
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

#-----------------------------初始条件------------------------------------------
E_ej  = 2e51                                           #爆发能量erg
M_ej  = 14*con.M_sun.to('g').value                     #爆发质量g
R_ej  = 1*con.pc.to('cm').value                        #爆发尺度cm
B     = 1.0e-6                                         #背景磁场G
lamda = 1.0e17                                         #最大MHD湍流波长
n_h   = 20.0                                           #H原子密度cm^-3
T     = 1.0e2                                          #温度
#t     = np.linspace(1,100000,10000)                    #时间yr
#-----------------------------参量设定------------------------------------------
radius= 37                                             #pc
pixel = 5120                                            #number
m_h   = con.m_p.to('g').value                          #H原子质量g
e     = con.e.esu.value                                #电子电荷C
u     = 1.3                                            #星际介质平均原子质量
s     = 1                                              #壳层速度分布
gamma = 5/3                                            #绝热指数
w_c   = 0.01                                            #外部幂率区域质量比重
n     = 0                                              #ejecta密度分布幂指数，Ia：n=7,II:n=9
#-----------------------------推导结果------------------------------------------
rho0  = n_h*u                                          #真实质量密度
r_c   = w_c*R_ej                                       #中心区半径
fn    = 3/4/ma.pi*(1-n/3)/(1-n/3*w_c**(3-n))           #质量守恒系数
f0    = fn*w_c**-n                                     #中心区密度
alpha = (3-n)/(5-n)*(w_c**(n-5)-n/5)/(w_c**(n-3)-n/3)*w_c**2
v_ej  = (E_ej/(M_ej*alpha*0.5))**0.5                   #最外层速度，cm/s
t     = (R_ej/v_ej*un.s).to('yr')                      #当前时间
rho_ch= M_ej/R_ej**3                                   #特征密度
R_ch  = (M_ej**(1.0/3.0)*(rho0*m_h)**(-1.0/3.0)*un.cm).to('pc').value

phi   = 0.55 
l     = 1.19
up    = 1 + (n-3)/3*(phi/l*fn)**0.5*((4*R_ej*un.cm).to('pc').value/R_ch)**1.5
down  = 1 + (n/3)*(phi/l*fn)**0.5*((4*R_ej*un.cm).to('pc').value/R_ch)**1.5
eta   = up/down
#==========================产生初始条件=========================================
r     = np.linspace(0,radius,pixel)
v     = np.linspace(0,radius,pixel)*0
rho   = np.linspace(0,radius,pixel)*0+rho0*m_h
p     = np.linspace(0,radius,pixel)*0

for i in range(len(r)):
    if r[i] <= 1: 
#        up    = 1 + (n-3)/3*(phi/l*fn)**0.5*(r[i]/R_ch)**1.5
#        down  = 1 + (n/3)*(phi/l*fn)**0.5*(r[i]/R_ch)**1.5
#        eta   = 1
        v[i]   = (r[i]*un.pc/t).to('cm/s').value
        if r[i] <= (r_c*un.cm).to('pc').value:
            rho[i] = rho_ch*fn*w_c**-n
        else:
            rho[i] = rho_ch*fn*(v[i]/v_ej)**-n
               
def fun_M(r):
    return 4*ma.pi*r**2*rho_ch*fn*(r/R_ej)**-n

def fun_E(r):
    return 4*ma.pi*r**2*rho_ch*fn*(r/R_ej)**-n*0.5*(r/t.to('s').value)**2
#==============================================================================
#fig, axes = plt.subplots(nrows=1, ncols=2)
#ax0, ax1 = axes.flat
#ax0.plot(r,v)
#ax0.set_xlim(0.0,1.1)
#ax0.set_xlabel('log(radius/pc)')
#ax0.set_ylabel('log(velocity/(km/s))')
#ax1.plot(r,np.log(rho))
#ax1.set_xlim(0.0,1.1)
#ax1.set_xlabel('log(radius/pc)')
#ax1.set_ylabel('log(density/(cm^-3))')