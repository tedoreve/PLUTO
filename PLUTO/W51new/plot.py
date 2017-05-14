import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp
import scipy.ndimage as nd

choose = 'single' #single or multiple or temp

wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

if choose == 'single':
    #D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
    D = pp.pload(30,w_dir=wdir)
#    print D.bx1.shape
    
    I = pp.Image()
    flux=D.rho*(D.bx1**2+D.bx2**2)**1.25*(D.vx1**2+D.vx2**2)**0
    print flux.shape
    flux= (flux-np.mean(flux))*1+np.mean(flux)*1
    flux=nd.gaussian_filter(flux,sigma=(14,14),order=0)
#    I.pldisplay(D, np.log10(flux),x1=D.x1, \
#                x2=D.x2,label1='l offset (pc)',label2='b offset (pc)',                                    \
#                title='Relative Radio Flux',                        
#                cbar=(True,'vertical'))
    imshow(np.log10(flux).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
    levels = np.arange(-7.1, -3.1, 0.5)
    contour(np.log10(flux).T, levels,
                 origin='lower',
                 linewidths=2,
                 extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
    cbar=colorbar()
    cbar.set_label('log(S)')
    xlabel('l offset (pc)')
    ylabel('b offset (pc)')
    title('Relative Radio Flux')
#    
#    T = pp.Tools()
#    newdims = 2*(20,)
#    Xmesh, Ymesh = meshgrid(D.x1.T,D.x2.T)
#    xcong = T.congrid(Xmesh,newdims,method='linear')
#    ycong = T.congrid(Ymesh,newdims,method='linear')
#    velxcong = T.congrid(D.bx1.T,newdims,method='linear')
#    velycong = T.congrid(D.bx2.T,newdims,method='linear')
#    gca().quiver(xcong, ycong, velxcong, velycong,color='w')
#    show()
    savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
#    close()

if choose == 'multiple':
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
        Xmesh, Ymesh = meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.vx1.T,newdims,method='linear')
        velycong = T.congrid(D.vx2.T,newdims,method='linear')
        gca().quiver(xcong, ycong, velxcong, velycong,color='w')
        savefig('MHD_Blast_'+str(i)+'.png') # Only to be saved as either .png or .jpg
        close()
#    show()
if choose == 'temp':
    D = pp.pload(30,w_dir=wdir)
    
    I = pp.Image()
    flux=D.prs*1.67e-7/D.rho/1000000/1.3806488e-23
    print flux.shape
#    flux= (flux-np.mean(flux))*5+np.mean(flux)*5.3
#    flux=nd.gaussian_filter(flux,sigma=(4,4),order=0)
    I.pldisplay(D, flux,x1=D.x1, \
                x2=D.x2,label1='l offset (pc)',label2='b offset (pc)',                                    \
                title='Temperature',                        
                cbar=(True,'vertical'))
#    savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
    show()
