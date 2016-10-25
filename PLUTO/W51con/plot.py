import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp
import scipy.ndimage as nd
wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

#D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
D = pp.pload(8,w_dir=wdir)
#print D.rho.shape[0]

#f=open('rho.txt','r+')
#temp=f.read()
#temp=

I = pp.Image()
flux=D.rho*(D.bx1**2+D.bx2**2)**2.0*(D.vx1**2+D.vx2**2)**0
print flux.shape
flux= (flux-np.mean(flux))*5+np.mean(flux)*5.1
flux=nd.gaussian_filter(flux,sigma=(15,15),order=0)
I.pldisplay(D, np.log(flux),x1=D.x1, \
            x2=D.x2,label1='x',label2='y',                                    \
<<<<<<< HEAD
            title=r'Flux [M=15 E=0.8 R=1 t=10 rho=50]',                         \
=======
            title=r'Flux [M=15 E=0.8 R=1 t=3 rho=20]',                         \
>>>>>>> 0c15addce092344029f364f8d9e72bba75d9063f
            cbar=(True,'vertical'))
#savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
show()
