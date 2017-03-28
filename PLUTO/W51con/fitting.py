# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 20:24:45 2017

@author: tedoreve
"""

import numpy as np  
import pylab as plt  
#import matplotlib.pyplot as plt  
from scipy.optimize import curve_fit  
from scipy import asarray as ar,exp  

def gaussian(x,*param):  
    result = param[0]*np.exp(-np.power(x - param[2], 2.) / (2 * np.power(param[4], 2.)))
    +param[1]*np.exp(-np.power(x - param[3], 2.) / (2 * np.power(param[5], 2.)))  
    
    return result
  
  
popt,pcov = curve_fit(gaussian,x,y,p0=[2,12000,0.01,2,130000,0.01])  
#print popt  
#print pcov  
  
plt.plot(x,y,'b+:',label='data')  
plt.plot(x,gaussian(x,*popt),'ro:',label='fit')  
plt.legend()  
plt.show()  
