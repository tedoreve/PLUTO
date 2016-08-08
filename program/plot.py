# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:03:37 2015

@author: lingzhi

"""

import matplotlib.pyplot as pl  
import numpy as np  
import subprocess  
import re
from mpl_toolkits.mplot3d import Axes3D

rho=np.loadtxt('rho.txt')
B=np.loadtxt('B.txt')

Eelec=np.loadtxt('Eelec.txt')
Felec=np.loadtxt('Felec.txt')
Ep=np.loadtxt('Ep.txt')
Fp=np.loadtxt('Fp.txt')
Ebk=np.loadtxt('Ebk.txt')
Fbk=np.loadtxt('Fbk.txt')



Eg=np.loadtxt('Eg.txt')

span=np.zeros(len(Eg))
span[99]=Eg[99]-Eg[98]
for i in range(len(Eg)-1):
    span[i]=Eg[i+1]-Eg[i]
F=np.zeros([100,100])
FF=np.zeros([100,100])
#for i in range(1000):
#    Fsyn=np.loadtxt('Fsyn'+str(i+1)+'.txt')
    #pl.plot(np.log10(sum(Eg[10:15])),np.log10(sum(Eg[10:15]*span[10:15]*Fsyn[10:15])))
#    F[i]=sum(Eg[10:15]*span[10:15]*Fsyn[10:15])
#pl.plot(rho,F)
for j in range(100):
	for i in range(100):
	    if i < 9:
		Fsyn=np.loadtxt('Fsyn'+str(j+1)+'           '+str(i+1)+'.txt')
	    elif 9 <= i < 99:
		Fsyn=np.loadtxt('Fsyn'+str(j+1)+'          '+str(i+1)+'.txt')
	    else: 
		Fsyn=np.loadtxt('Fsyn'+str(j+1)+'         '+str(i+1)+'.txt')
	    #pl.plot(np.log10(sum(Eg[10:15])),np.log10(sum(Eg[10:15]*span[10:15]*Fsyn[10:15])))
	    F[j,i]=sum(Eg[10:15]*span[10:15]*Fsyn[10:15])
            FF[j,i]=sum(Eg[55:60]*span[55:60]*Fsyn[55:60])
fig=pl.figure()
ax = fig.add_subplot(111,projection='3d')
[x,y]=np.meshgrid(rho,B)
ax.scatter(x,y,F,c='w',marker='o')
ax.scatter(x,y,FF,c='r',marker='o')

Fbrem=np.loadtxt('Fbrem.txt')
Fpp=np.loadtxt('Fpp.txt')
Fic=np.loadtxt('Fic.txt')

#x1=np.loadtxt('x1.txt')
#x2=np.loadtxt('x2.txt')
#x3=np.loadtxt('x3.txt')


#fig1=pl.figure('fig1')
#pl.plot(np.log10(Eelec),np.log10(Eelec**2*Felec),label='electron')
#pl.plot(np.log10(Ep),np.log10(Ep**2*Fp),label='proton')
#pl.plot(np.log10(Ebk/10e9),np.log10(Ebk**2*Fbk*10e9),label='background')
#pl.ylim(-30.-7)
#pl.legend(loc='best')
#pl.savefig('/home/lingzhi/Downloads/program/1.png')
#print span
#fig2=pl.figure('fig2')
 

#pl.plot(np.log10(Eg),np.log10(Eg*span*Fbrem),label='brem')
#pl.plot(np.log10(Eg),np.log10(Eg*span*Fpp),label='pp')
#pl.plot(np.log10(Eg),np.log10(Eg*span*Fic),label='ic')
#pl.plot(np.log10(Eg),np.log10(Eg**2*Fic+Eg**2*Fpp+Eg**2*Fbrem+Eg**2*Fsyn),label='total')
#pl.ylim(-14,-7)
#pl.xlabel('$log(E),(GeV)$')
#pl.ylabel('$log(E^2\cdot F),(GeV\cdot cm^{-3}\cdot s^{-1})$')
#pl.legend(loc='best')
#pl.savefig('/home/lingzhi/Downloads/program/2.png')

#for i in range(len(x1)-1):
#delta=x1[2]-x1[1]
#pl.plot(np.log10(x1/1.0e6),np.log10(delta*x2*x1/1.0e6),label='syn')
#pl.plot(np.log10(x1/1.0e6),np.log10(delta*x3*x1/1.0e6),label='brem')

#pl.ylim(-20,0)
#pl.xlabel('$log(E),(GeV)$')
#pl.ylabel('$log(E^2\cdot F),(GeV\cdot cm^{-3}\cdot s^{-1})$')
#pl.legend(loc='best')
pl.show()
#pl.savefig('/home/lingzhi/Downloads/program/3.png')


