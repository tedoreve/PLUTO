import os
import sys
from numpy import *
from matplotlib.pyplot import *
import pyPLUTO as pp
from mpl_toolkits.mplot3d import Axes3D

wdir = '/home/lingzhi/Downloads/Github/PLUTO/PLUTO/BSNR/'
nlinf = pp.nlast_info(w_dir=wdir)

#D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
D = pp.pload(1,w_dir=wdir)
print D.rho.shape

fig=figure()
ax = fig.add_subplot(111,projection='3d')
[x,y,z]=np.meshgrid(D.x1,D.x2,D.x3)
ax.scatter(x, y, z, c='r', marker='o')

#f=open('rho.txt','r+')
#temp=f.read()
#temp=

#I = pp.Image()
#I.pldisplay(D, log10(D.rho),x1=D.x1,x2=D.x2,x3=60,label1='x',label2='y',title=r'Density $\rho$ [MHD_Blast test]',cbar=(True,'vertical'))
#savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
show()
