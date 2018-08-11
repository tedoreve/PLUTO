
import numpy as np
import matplotlib.pyplot as plot
import pyPLUTO as pp
import zmf as z
import astropy.constants as con
import astropy.units as un

def cf(temperature,plot=False):
    '''
    calculate cooling function lambda(T)
    temperature(K), return(erg cm^-3 s^-1)
    array input is accessible
    '''
    cooling = np.loadtxt('cooling.dat')
    if plot:
        plt.loglog(cooling[:,0],cooling[:,1])
        plt.loglog(temperature,np.interp(temperature, cooling[:,0],cooling[:,1]),'o')
    
    return np.interp(temperature, cooling[:,0],cooling[:,1])

wdir = './'
nlinf = pp.nlast_info(w_dir=wdir)

#D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
D = pp.pload(2,w_dir=wdir)
T= (D.prs*un.Ba/(D.rho*un.cm**-3)/con.k_B*1.67e-6).to('K')
print cf(T.value).shape
#f=open('rho.txt','r+')
#temp=f.read()
#temp=
#imshow(T.value)

I = pp.Image()
I.pldisplay(D, cf(T.value),x1=D.x1,x2=D.x2,label1='x',label2='y',title=r'Density $\rho$ [MHD_Blast test]',cbar=(True,'vertical'))
#savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
plt.show()
