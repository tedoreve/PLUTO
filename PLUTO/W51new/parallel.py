# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 22:34:33 2017

@author: tedoreve
"""

import numpy as np
from multiprocessing import Pool
import time

#def toff(f):
#    def wrapper(*args):
#        start = time.time()
#        f(*args)
#        end   = time.time()
#        print(end-start)
#    return wrapper
#
#@toff
def f(i,j,k):
    return i+j+k


if __name__=='__main__':
    width = 512
#    x   = np.zeros([width,width,width])
#    p = Pool(4)
    start = time.time()
    x = np.fromfunction(f,(width,width,width))
#    for i in range(width):
##        print(i)
#        for j in range(width):
#            for k in range(width):
#                x[i,j,k] = i + j + k
#            for l in range(2*width):
#                print(p.apply_async(judge, args=(i,j,l,x)).get())

#    p.close()
#    p.join()
    end = time.time()
    print('Total time.', end - start)
    print(x)


#@toff
#def magnetism(width,widthi,widthj):
#width = 128
#bx1 = np.zeros([width,width])
#bx2 = np.zeros([width,width])
#x   = np.zeros([width,width])
#
#p = Pool(64)
#
#for i in range(width):
#    for j in range(width):
#        for l in range(2*width):
#            x[i,j] = p.apply_async(judge, args=(i,j,l)).get()
#p.close()
#p.join()                 
#
#
#bx1  = np.rot90(x)
#bx2  = np.rot90(x)
#
#bx1 = bx1.T/500000
#bx2 = bx2.T/500000
#
#bx1 = np.reshape(bx1,width**2,1)
#bx2 = np.reshape(bx2,width**2,1)
#    return bx1 , bx2
#

#bx1 , bx2 = magnetism(width,width,width)
