import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')   # generate postscript output by default
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import pyPLUTO as pp
import scipy.ndimage as nd

def single(ty,t,E,rho,sigma,wdir):
    #D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
    D = pp.pload(t/1000-1,w_dir=wdir)   
    flux=D.rho*(D.bx1**2+D.bx2**2)**1.25*(D.vx1**2+D.vx2**2)**0
    flux= (flux-np.mean(flux))*1+np.mean(flux)*1
    flux=nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
    if ty == 'flux':
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        neg = ax.imshow(np.log10(flux).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        levels = np.arange(-5.5, -4, 0.5)
        plt.contour(np.log10(flux).T, levels,
                     origin='lower',
                     linewidths=2,
                     extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        cbar = fig.colorbar(neg,ax=ax)
        cbar.set_label(r'log(S)')
        ax.set_xlabel('l offset (pc)')
        ax.set_ylabel('b offset (pc)')
        ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
        fig.savefig('t='+str(t)+' density='+str(rho)+' E='+str(E)+'.eps') # Only to be saved as either .png or .jpg

    else:
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        neg = ax.imshow(np.log10(D.rho).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        cbar = fig.colorbar(neg,ax=ax)
        cbar.set_label(r'log($\mathregular{\rho/cm^{-3}}$)')
        ax.set_xlabel('l offset (pc)')
        ax.set_ylabel('b offset (pc)')
        ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx1.T,newdims,method='linear')
        velycong = T.congrid(D.bx2.T,newdims,method='linear')
        plt.gca().quiver(xcong, ycong, velxcong, velycong,color='w')
    #    plt.show()
        fig.savefig('rho-t='+str(t)+' density='+str(rho)+' E='+str(E)+'.eps') # Only to be saved as either .png or .jpg
    #    close()
        

def multiple(wdir):
    #D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
    for i in range(30):
        D = pp.pload(i,w_dir=wdir)
        #print D.rho.shape[0]
        
        I = pp.Image()
        I.pldisplay(D, np.log(D.rho),x1=D.x1, \
                    x2=D.x2,label1='l offset (pc)',label2='b offset (pc)',                                    \
                    title=r'$Density (ln(cm^{-3})) - Velocity (km/s)$ '+str(i*1000)+' years',    
                    cbar=(True,'vertical'))
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.vx1.T,newdims,method='linear')
        velycong = T.congrid(D.vx2.T,newdims,method='linear')
        plt.gca().quiver(xcong, ycong, velxcong, velycong,color='w')
        plt.savefig('20_'+str(i)+'.png') # Only to be saved as either .png or .jpg
        plt.close()
#    show()
def temp(wdir):
    D = pp.pload(30,w_dir=wdir)
    
    I = pp.Image()
    flux=D.prs*1.67e-7/D.rho/1000000/1.3806488e-23
#    print flux.shape
#    flux= (flux-np.mean(flux))*5+np.mean(flux)*5.3
#    flux=nd.gaussian_filter(flux,sigma=(4,4),order=0)
    I.pldisplay(D, flux,x1=D.x1, \
                x2=D.x2,label1='l offset (pc)',label2='b offset (pc)',                                    \
                title='Temperature',                        
                cbar=(True,'vertical'))
#    savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
    plt.show()
#==================================main========================================  
if __name__=='__main__':   
    #font = {'family' : 'serif',  
    #        'weight' : 'normal',  
    #        'size'   : 12,  
    #        }  
    
    choose = 'single' #single or multiple or temp
    ty     = 'flux'    
    t      = 18000   
    E      = 1.3
    rho    = 0.3
    sigma  = 14
    
    wdir = './'
    nlinf = pp.nlast_info(w_dir=wdir)
    single(ty,t,E,rho,sigma,wdir)
