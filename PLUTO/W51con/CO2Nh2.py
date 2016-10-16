import numpy as np
import matplotlib.pyplot as plt
#matplotlib.use('Agg')

from matplotlib.pyplot import plot,savefig
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from astropy.io import fits
from astropy import wcs
import warnings
from astropy.convolution import convolve,Gaussian2DKernel,Tophat2DKernel
from astropy.modeling.models import Gaussian2D
from astropy import units as un
from astropy import constants as con
import astropy.visualization as v
import astropy.tests.zmf as z

hdulist = fits.open('grs-49-cube.fits')
hdulist.info()                         #显示fits信息 
pchead = hdulist[0].header#读取fits文件头
pcdata = hdulist[0].data  #读取fits数据   
hdulist.close()
w=wcs.WCS(pchead)

l1=49.70
l2=48.70
b1=-0.95
b2=0.05

x1,y1,a=w.wcs_world2pix(l1,b1,0,0)         #用读取的坐标数据将经纬度转变成pixel数
x2,y2,a=w.wcs_world2pix(l2,b2,0,0)
img=(pcdata[354,:,:])*2e20/(z.distance2diameter(4.3,60)*un.pc.to('cm'))
img=np.abs(np.nan_to_num(img))
plt.imshow(img[y1:y2,x1:x2],origin='lower',interpolation='nearest',      \
               extent=[l1,l2,b1,b2])
#               
xx=img[y1:y2,x1:x2]

pchead['NAXIS']=2
pchead['NAXIS1']=162
pchead['NAXIS2']=162
pchead['CRVAL1']=49.2
pchead['CRPIX1']=81
pchead['CRVAL2']=-0.45
pchead['CRPIX2']=81
pchead.remove('NAXIS3')
#pchead.remove('NAXIS4')
pchead.remove('CTYPE3')
pchead.remove('CRVAL3')
pchead.remove('CRPIX3')
pchead.remove('CDELT3')
pchead.remove('CROTA3')
#pchead.remove('CTYPE4')
#pchead.remove('CRVAL4')
#pchead.remove('CRPIX4')
#pchead.remove('CDELT4')
#pchead.remove('CROTA4')
#pchead.remove('DATAMIN')
#pchead.remove('DATAMAX')
#pchead.remove('MINCOL')
#pchead.remove('MAXCOL')
#pchead.remove('MINROW')
#pchead.remove('MAXROW')
#pchead.remove('MINFIL')
#pchead.remove('MAXFIL')


prihdu=fits.PrimaryHDU(header=pchead,data=xx)
hdu=fits.HDUList([prihdu])
hdu.writeto('W51C.fits')

