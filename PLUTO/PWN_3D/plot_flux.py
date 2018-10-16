# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 12:44:32 2018

@author: tedoreve
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyPLUTO as pp
import numpy as np
import scipy.ndimage as nd
from mpl_toolkits.axes_grid1 import AxesGrid

#==============================================================================

def plot(D,i,flux_vol,xlabel,ylabel,name,sigma):
    
    flux = flux_vol.sum(axis=i)
    flux = nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
    
    fig = plt.figure(figsize=(6.5,6.0))
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

    neg = ax[0].imshow(flux.T,origin='lower',extent=[D.x2[0],D.x2[-1],D.x3[0],D.x3[-1]])

    
    ax.cbar_axes[0].colorbar(neg)
    
    
#    fig.subplots_adjust(top=0.99,bottom=0.11,left=0.11,right=0.99)
    
    ax.axes_llc.set_xlabel(xlabel,fontsize=20)
    ax.axes_llc.set_ylabel(ylabel,fontsize=20)       
    fig.savefig(name)
    

#flux_vol = (flux-np.mean(flux))*1+np.mean(flux)*1
#==============================================================================
if __name__ == '__main__':
    t      = 4 
    sigma  = 1.0
    vindex = 0.0
    
    wdir = './'
    nlinf = pp.nlast_info(w_dir=wdir)
    
    #==============================================================================
    
    D = pp.pload(t,w_dir=wdir)   
    v_mean = np.mean(D.vx1**2+D.vx2**2+D.vx3**2)
    print(np.max(D.rho),np.min(D.rho),np.max((D.bx1**2+D.bx2**2+D.bx3**2)**0.5),
          np.min((D.bx1**2+D.bx2**2+D.bx3**2)**0.5))
    #print(np.sign(D.vx1**2+D.vx2**2+D.vx3**2-v_mean))
    for i in range(3):
    
        
        if i == 0:
            xlabel = 'y (pc)'
            ylabel = 'z (pc)'
            name   = 't'+str(t)+'_yz.eps'
    #        flux_vol = np.log10(D.prs/D.rho)*(D.bx2**2+D.bx3**2)**0.0*1e5**(-(np.sign(D.vx1**2+D.vx2**2+D.vx3**2-v_mean)+1)*(D.vx1**2+D.vx2**2+D.vx3**2)/2)
                    
            flux_vol = D.prs/D.rho
    
        if i == 1:
            xlabel = 'x (pc)'
            ylabel = 'z (pc)'
            name   = 't'+str(t)+'_xz.eps'
            #flux_vol = D.rho*((D.bx1/2**0.5)**2+(D.bx2/2**0.5)**2+D.bx3**2)**2.0*(D.vx1**2+D.vx2**2+D.vx3**2)**vindex
            
            flux_vol = D.prs/D.rho
        if i == 2:
            xlabel = 'x (pc)'
            ylabel = 'y (pc)'
            name   = 't'+str(t)+'_xy.eps'
            flux_vol = D.rho*(D.bx1**2+D.bx2**2)**2.0*(D.vx1**2+D.vx2**2+D.vx3**2)**vindex
            
        plot(i,flux_vol,xlabel,ylabel,name,sigma)

