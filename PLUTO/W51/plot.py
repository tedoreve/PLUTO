import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp

wdir = '/home/lingzhi/Downloads/Github/PLUTO/PLUTO/W51/'
nlinf = pp.nlast_info(w_dir=wdir)

D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
#D = pp.pload(0,w_dir=wdir,datatype='float')
#print D.rho.shape[0]

#f=open('rho.txt','r+')
#temp=f.read()
#temp=

I = pp.Image()
I.pldisplay(D, D.rho,x1=D.x1,x2=D.x2,label1='x',label2='y',title=r'Density $\rho$ [MHD_Blast test]',cbar=(True,'vertical'))
#savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
show()
