
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
E_ej  = 2e51                                           #爆发能量erg
M_ej  = 14*con.M_sun.to('g').value                     #爆发质量g
R_ej  = 1*con.pc.to('cm').value                        #爆发尺度cm
B     = 1.0e-6                                         #背景磁场G
lamda = 1.0e17                                         #最大MHD湍流波长
n_h   = 20.0                                           #H原子密度cm^-3
T     = 1.0e2                                          #温度
t     = np.linspace(1,100000,10000)                   #时间yr
#-----------------------------参量设定------------------------------------------
m_h   = con.m_p.to('g').value                          #H原子质量g
e     = con.e.esu.value                                #电子电荷C
u     = 1.3                                            #星际介质平均原子质量
s     = 1                                              #壳层速度分布
gamma = 5/3                                            #绝热指数
eta   = 3/7                                            #外部幂率区域质量比重
n     = 2                                              #ejecta密度分布幂指数，Ia：n=7,II:n=9
#-----------------------------推导结果------------------------------------------
c_s     = 9.2e5*np.sqrt(gamma*T/1.0e4)                 #当地声速cm/s
n_ISM   = n_h+0.3*2*n_h                                #星际介质密度cm^-3
rho_ISM = n_ISM*m_h                                    #星际介质密度g/cm^-3
rho_ISM = rho_ISM                                      #介质密度分布
rho_0   = rho_ISM                                      #R_ej处密度
r_c     = ((1-eta*n/3)/(1-eta))**(1/(n-3))*R_ej        #内ejecta均匀介质半径
rho_c   = (1-eta)*M_ej/(4/3*ma.pi*r_c**3)              #内ejecta均匀介质密度
rho_0   = rho_c*(R_ej/r_c)**(-n)                       #ejecta外边界密度
E       = E_ej/(4/3*ma.pi*R_ej**3)                     #能量密度
v0      = ma.sqrt(E_ej)/ma.sqrt(2*ma.pi*rho_c*r_c**5/5/R_ej**2+2*ma.pi*rho_0*      \
          R_ej**3*(1-(R_ej/r_c)**(n-5))/(5-n))         #ejecta外边界速度
P       = n_ISM*con.k_B.to('erg/K').value*T            #压强
t_Sed   = (3*M_ej/4/np.pi/u/m_h/n_h/v0**3)**(1/3)      #s,可用t_Sed*un.s.to('yr')
t_Sed   = t_Sed*un.s.to('yr')
t_Rad   = 2.7e4*(E_ej/1.0e51)**0.24*n_h**(-0.52)       #yr     

C       = 2*ma.pi*r_c**n*rho_c*R_ej**(-2*s)

index   = 2*s+3-n
#E_out   = C/index*(r_c**index-R_ej**index)
#E_inner = 2*ma.pi*rho_c*v0**2*R_ej**(-2*s)*r_c**(2*s+3)/(2*s+3)
#print(E_out,E_inner)

part1   = C/r_c**n*r_c**(2*s+3)/(2*s+3)
part2   = C*(R_ej**index-r_c**index)/index
v0      = E_ej**0.5*(part1+part2)**(-0.5)

#------------------------------演化结果-----------------------------------------
f       = 10                                           #粒子平均自由程/回旋半径
Rj      = 1                                            #激波相对磁场方向
v=np.zeros(len(t))
R=np.zeros(len(t))
E_e=np.zeros(len(t))
E_p=np.zeros(len(t))
E_max1=np.zeros(len(t))
E_max2=np.zeros(len(t))
E_max3=np.zeros(len(t))

#下面计算最大加速能量
def fun(t,v,B,f,Rj):
    return 100*B*(v/1.0e8)**2/f/Rj

for i in range(len(t)):
    if t[i] < t_Sed:
        v[i]      = v0                                              #t<t_Sed,yr
        R[i]      = v0*t[i]*3600*24*365                             #t<t_Sed,yr
        E_max1[i] = sc.integrate.quad(fun, 0, t[i]*3600*24*365, args=(v[i],B,f,Rj))[0]
    if t_Sed < t[i] < t_Rad:
        v[i]      = v0*(t_Sed/t[i])**0.6                            #t_Sed<t<t_Rad,yr
        R[i]      = 2.5*v0*3600*24*365*t_Sed*((t[i]/t_Sed)**0.4-0.6)#t_Sed<t<t_Rad,yr
        E_max1[i] = sc.integrate.quad(fun, 0, t_Sed*3600*24*365, args=(v0,B,f,Rj))[0]+\
                    sc.integrate.quad(fun, t_Sed*3600*24*365, t[i]*3600*24*365, args=(v[i],B,f,Rj))[0]
    if t[i] > t_Rad:
        v[i]      = v0*(t_Sed/t_Rad)**0.6*(t_Rad/t[i])**0.69        #t_Rad<t,yr
        R[i]      = 2.5*v0*3600*24*365*t_Sed*(1.29*(t_Rad/t_Sed)**0.4*((t[i]/t_Rad)**0.31-0.225)-0.6)   #t_Rad<t,yr
        E_max1[i] = sc.integrate.quad(fun, 0, t_Sed*3600*24*365, args=(v0,B,f,Rj))[0]+\
                    sc.integrate.quad(fun, t_Sed*3600*24*365, t_Rad*3600*24*365, args=(v[i],B,f,Rj))[0]+\
                    sc.integrate.quad(fun, t_Rad*3600*24*365, t[i]*3600*24*365, args=(v[i],B,f,Rj))[0]
             
    E_max2[i]  = 2e5*(f*Rj*B)**(-0.5)*(v[i]/1.0e8)
    E_max3[i]  = e*B*lamda/f*un.erg.to('MeV')
    E_e[i]     = min(E_max1[i],E_max2[i],E_max3[i])    
    E_p[i]     = min(E_max1[i],E_max3[i])    
#------------------------------激波相关----------------------------------------
M       = v/c_s                                        #马赫数
alpha   = (gamma+1)*M**2/((gamma-1)*M**2+2)            #压缩比
a       = (2+alpha)/(alpha-1)                          #粒子谱指数
E_gain  = 100*B*(v/1.0e8)**2/f/Rj                      #粒子在激波中的能量获得率，MeV/s 
#plt.plot(t,E)
#plt.plot(np.log10(t),np.log10(E))
#plt.plot(np.log10(t),np.log10(v))
fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flat
ax0.plot(np.log10(t),np.log10(v))
ax0.set_xlabel('log(time/year)')
ax0.set_ylabel('log(velocity/(cm/s))')
ax1.plot(np.log10(t),np.log10(R*un.cm.to('pc')))
ax1.set_xlabel('log(time/year)')
ax1.set_ylabel('log(radius/pc)')
ax2.plot(np.log10(t),np.log10(E_e))
ax2.set_xlabel('log(time/year)')
ax2.set_ylabel('log(E$_e$/MeV)')
ax3.plot(np.log10(t),np.log10(E_p))
ax3.set_xlabel('log(time/year)')
ax3.set_ylabel('log(E$_p$/MeV)')


    