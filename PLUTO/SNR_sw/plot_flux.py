# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 12:44:32 2018

@author: tedoreve
"""


import matplotlib.pyplot as plt
import pyPLUTO as pp
import scipy.ndimage as nd
from mpl_toolkits.axes_grid1 import AxesGrid

#==============================================================================

ty     = 'flux'    
t      = 1 
E      = 1.3
rho    = 1.0
sigma  = 1

wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

#==============================================================================

D = pp.pload(t,w_dir=wdir)   


flux_vol = D.rho*(D.bx1**2+D.bx2**2+D.bx3**2)**1.25*(D.vx1**2+D.vx2**2+D.vx3**2)**0
#flux_vol = (flux-np.mean(flux))*1+np.mean(flux)*1

for i in range(3):
    flux = flux_vol.sum(axis=i)
    flux = nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
    
    fig = plt.figure(figsize=(6.3,5.5))
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

    neg = ax[0].imshow(flux.T*10000,origin='lower',extent=[D.x2[0],D.x2[-1],D.x3[0],D.x3[-1]])

    
    ax.cbar_axes[0].colorbar(neg)
    
    ax.cbar_axes[0].toggle_label(True)
    ax.cbar_axes[0].axis[ax.cbar_axes[0].orientation].set_label('relative flux density')
    
    fig.subplots_adjust(top=0.99,bottom=0.08,left=0.11,right=0.91)
    
    if i == 0:
        ax.axes_llc.set_xlabel('y (pc)')
        ax.axes_llc.set_ylabel('z (pc)')       
        fig.savefig('t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_yz.eps')
        
    if i == 1:
        ax.axes_llc.set_xlabel('x (pc)')
        ax.axes_llc.set_ylabel('z (pc)')       
        fig.savefig('t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_xz.eps')
        
    if i == 2:
        ax.axes_llc.set_xlabel('x (pc)')
        ax.axes_llc.set_ylabel('y (pc)')       
        fig.savefig('t'+str(t)+'_density'+str(int(rho))+'_E'+str(int(E))+'_xy.eps')
