
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:12:01 2016

@author: tedoreve
"""

from astropy import constants as con
from astropy import units as un
import math as ma
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

#-----------------------------初始条件------------------------------------------
E_ej  = 1#e51                                           #爆发能量erg
M_ej  = 14#*con.M_sun.to('g').value                     #爆发质量g
R_ej  = 1#*con.pc.to('cm').value                        #爆发尺度cm
B     = 1.0e-6                                         #背景磁场G
lamda = 1.0e17                                         #最大MHD湍流波长
n_h   = 0.4                                           #H原子密度cm^-3
T     = 1.0e2                                          #温度
#-----------------------------参量设定------------------------------------------
m_h   = con.m_p.to('g').value                          #H原子质量g
e     = con.e.esu.value                                #电子电荷C
u     = 1.3                                            #星际介质平均原子质量
s     = 0                                              #壳层速度分布
gamma = 5/3                                            #绝热指数
eta   = 3/7                                            #外部幂率区域质量比重
n     = 9                                              #ejecta密度分布幂指数，Ia：n=7,II:n=9
#-----------------------------推导结果------------------------------------------
c_s     = 9.2e5*np.sqrt(gamma*T/1.0e4)                 #当地声速cm/s

#P       = n_ISM*con.k_B.to('erg/K').value*T            #压强
   
#-------------------------------类高斯函数-------------------------------------
#构造一个类高斯函数
def g1(x,u,o1,o2):
    return np.exp(-0.5*((x-u)/((np.arctan(x-u)+np.pi)/2/np.pi*(o2-o1)+o1))**2)
def g2(x,u,o1,o2):
    return np.exp(-0.5*((x-u)/((o2-o1)/(1+e**(x-u))+o1))**2)
def g3(x,u,o1,o2):
    return np.exp(-0.5*((x-u)/o2)**2)
x=np.linspace(0,R_ej*2,1000)
o1=R_ej/10
o2=R_ej/10
u=R_ej 
plt.plot(x,g2(x,u,o1,o2),label='o1=1,o2=1')  
#plt.ylim(-0.01,0.01)
plt.legend()
def fun_M(x,u,o1,o2):
    return 4*ma.pi*x**2*np.exp(-0.5*((x-u)/((o2-o1)/(1+e**(x-u))+o1))**2)

def fun_E(x,u,o1,o2):
    return 4*ma.pi*x**2*np.exp(-0.5*((x-u)/((o2-o1)/(1+e**(x-u))+o1))**2)*0.5*np.exp(-0.5*((x-u)/((o2-o1)/(1+e**(x-u))+o1))**2)**2

print(sc.integrate.quad(fun_E, 0, R_ej, args=(u,o1,o2))[0])
#------------------------------演化结果-----------------------------------------

#fig, axes = plt.subplots(nrows=2, ncols=2)
#ax0, ax1, ax2, ax3 = axes.flat
#ax0.plot(np.log10(t),np.log10(v))
#ax0.set_xlabel('log(time/year)')
#ax0.set_ylabel('log(velocity/(cm/s))')
#ax1.plot(np.log10(t),np.log10(R*un.cm.to('pc')))
#ax1.set_xlabel('log(time/year)')
#ax1.set_ylabel('log(radius/pc)')
#ax2.plot(np.log10(t),np.log10(E_e))
#ax2.set_xlabel('log(time/year)')
#ax2.set_ylabel('log(E$_e$/MeV)')
#ax3.plot(np.log10(t),np.log10(E_p))
#ax3.set_xlabel('log(time/year)')
#ax3.set_ylabel('log(E$_p$/MeV)')


    