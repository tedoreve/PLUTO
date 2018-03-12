# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 16:43:51 2018

@author: tedoreve
"""
import numpy as np
import matplotlib.pyplot as plt
import pyPLUTO as pp
from mpl_toolkits.axes_grid1 import AxesGrid

#==============================================================================
    
t      = 6 
E      = 1.3
rho    = 1.0
sigma  = 1
index  = 64

wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

#==============================================================================

D = pp.pload(t,w_dir=wdir)   


for i in range(3):
    
    fig = plt.figure(figsize=(6.2,5.5))
    ax = AxesGrid(fig, 111,  # similar to subplot(122)
            nrows_ncols=(1, 1),
            axes_pad=0.0,
            label_mode="1",
            share_all=True,
            cbar_location="right",
            cbar_mode="edge",
            cbar_size="2%",
            cbar_pad="0%",
            )
   
    fig.subplots_adjust(top=0.99,bottom=0.11,left=0.11,right=0.95)
    
    if i == 0:
        neg = ax[0].imshow(np.log10(D.rho[index,:,:]).T,origin='lower',extent=[D.x2[0],D.x2[-1],D.x3[0],D.x3[-1]])
        ax.cbar_axes[0].colorbar(neg)        
        
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x2.T,D.x3.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx2[index,:,:].T,newdims,method='linear')
        velycong = T.congrid(D.bx3[index,:,:].T,newdims,method='linear')
        ax[0].quiver(xcong, ycong, velxcong, velycong,color='w')
            
        ax.axes_llc.set_xlabel('y (pc)',fontsize=20)
        ax.axes_llc.set_ylabel('z (pc)',fontsize=20)       
        fig.savefig('rho_t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_yz.eps')
        
    if i == 1:
        neg = ax[0].imshow(np.log10(D.rho[:,index,:]).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x3[0],D.x3[-1]])
        ax.cbar_axes[0].colorbar(neg)        
        
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x3.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx1[:,index,:].T,newdims,method='linear')
        velycong = T.congrid(D.bx3[:,index,:].T,newdims,method='linear')
        ax[0].quiver(xcong, ycong, velxcong, velycong,color='w')
        
        ax.axes_llc.set_xlabel('x (pc)',fontsize=20)
        ax.axes_llc.set_ylabel('z (pc)',fontsize=20)       
        fig.savefig('rho_t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_xz.eps')
        
    if i == 2:
        neg = ax[0].imshow(np.log10(D.rho[:,:,index]).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        ax.cbar_axes[0].colorbar(neg)        
        
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx1[:,:,index].T,newdims,method='linear')
        velycong = T.congrid(D.bx2[:,:,index].T,newdims,method='linear')
        ax[0].quiver(xcong, ycong, velxcong, velycong,color='w')
        
        ax.axes_llc.set_xlabel('x (pc)',fontsize=20)
        ax.axes_llc.set_ylabel('y (pc)',fontsize=20)       
        fig.savefig('rho_t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_xy.eps')

