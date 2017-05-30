# -*- coding: utf-8 -*-
"""
Created on Sat May 27 09:16:18 2017

@author: tedoreve
"""

import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt

width = 128

def f(i,j,k):
    return i+j
    
def g(i,j):
    return i+j
    
#V = np.fromfunction(f,(width,width,width))/500000
#V = np.rot90(V,k=2).T
#S = np.fromfunction(g,(width,width))/500000


#    
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(V),
#                            plane_orientation='x_axes',
#                            slice_index=10,
#                        )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(V),
#                            plane_orientation='y_axes',
#                            slice_index=10,
#                        )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(V),
#                            plane_orientation='z_axes',
#                            slice_index=10,
#                        )
#
#mlab.outline()
index=2.4
rho_constant=0.13
pcdata = 1/np.random.power(index,[width,width,width])*rho_constant
#    pcdata = np.zeros([width,width,width])+rho_constant
for k in range(width):
    for j in range(width):
        for i in range(width):
            if (i-width/2)**2+(j-width/2)**2+(k-20)**2 < (width/8)**2:
                pcdata[k,j,i] = 20
    rho    = np.reshape(pcdata,width**3,1)
mlab.pipeline.volume(mlab.pipeline.scalar_field(pcdata))

#plt.imshow(S)