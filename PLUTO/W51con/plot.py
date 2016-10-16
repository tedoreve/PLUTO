import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

#D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
D = pp.pload(10,w_dir=wdir)
#print D.rho.shape[0]

#f=open('rho.txt','r+')
#temp=f.read()
#temp=

I = pp.Image()
flux=D.rho*(D.bx1**2+D.bx2**2)**2.0*(D.vx1**2+D.vx2**2)**0
print np.max(flux),np.min(flux)
result= (flux-np.mean(flux))*5+np.mean(flux)*6
I.pldisplay(D, np.log(result),x1=D.x1, \
            x2=D.x2,label1='x',label2='y',                                    \
            title=r'Flux [M=15 E=0.8 R=1 t=3 rho=100]',                         \
            cbar=(True,'vertical'))
#savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
show()
