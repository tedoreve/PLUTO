import numpy as np
import matplotlib as mpl
#mpl.use('Agg')   # generate postscript output by default
import matplotlib.pyplot as plt
#import matplotlib.pylab as plb
import pyPLUTO as pp
import scipy.ndimage as nd
#from matplotlib.ticker import ScalarFormatter
from mayavi import mlab
from matplotlib.colors import Normalize
from tvtk.util.ctf import ColorTransferFunction
from tvtk.util.ctf import PiecewiseFunction
from mpl_toolkits.axes_grid1 import AxesGrid


#class MidpointNormalize(Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        Normalize.__init__(self, vmin, vmax, clip)
#
#    def __call__(self, value, clip=None):
#        # I'm ignoring masked values and all kinds of edge cases to make a
#        # simple example...
#        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y))
#norm = MidpointNormalize(midpoint=0)
norm = Normalize(vmin=-3.5, vmax=-1.9)  

def single(ty,t,E,rho,sigma,wdir):
    #D = pp.pload(nlinf['nlast'],w_dir=wdir) # Loading the data into a pload object D
    D = pp.pload(t/1000-1,w_dir=wdir)   
    flux=D.rho*(D.bx1**2+D.bx2**2)**1.25*(D.vx1**2+D.vx2**2)**0
    flux= (flux-np.mean(flux))*1+np.mean(flux)*1
    flux=nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
    if ty == 'flux':
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        neg = ax.imshow(np.log10(flux).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        levels = np.arange(-5.5, -4, 0.5)
        plt.contour(np.log10(flux).T, levels,
                     origin='lower',
                     linewidths=2,
                     extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        cbar = fig.colorbar(neg,ax=ax)
        cbar.set_label(r'log(S)')
        ax.set_xlabel('l offset (pc)')
        ax.set_ylabel('b offset (pc)')
        ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
        fig.savefig('t='+str(t)+' density='+str(rho)+' E='+str(E)+'.eps') # Only to be saved as either .png or .jpg

    else:
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        neg = ax.imshow(np.log10(D.rho).T,origin='lower',extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
        cbar = fig.colorbar(neg,ax=ax)
        cbar.set_label(r'log($\mathregular{\rho/cm^{-3}}$)')
        ax.set_xlabel('l offset (pc)')
        ax.set_ylabel('b offset (pc)')
        ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx1.T,newdims,method='linear')
        velycong = T.congrid(D.bx2.T,newdims,method='linear')
        plt.gca().quiver(xcong, ycong, velxcong, velycong,color='w')
    #    plt.show()
        fig.savefig('rho-t='+str(t)+' density='+str(rho)+' E='+str(E)+'.eps') # Only to be saved as either .png or .jpg
    #    close()

def analyze(ty,t,E,rho,sigma,wdir):      
    D = pp.pload(t/1000-1,w_dir=wdir)   
    flux   = D.rho*(D.bx1**2+D.bx2**2)**1.25*(D.vx1**2+D.vx2**2)**0
    flux   = nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
    mid    = np.rot90(flux) 
    mid    = mid*np.eye(mid.shape[0])
    mid    = mid.sum(axis=0)
    
    
    fig = plt.figure(figsize=(7,6))
    ax = fig.add_subplot(111)
    ax.plot(D.x1*np.sqrt(2),mid)
    ax.set_xlabel('offset (pc)')
    ax.set_ylabel('S')
    ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
    ax.set_xlim(D.x1.min()*np.sqrt(2),D.x1.max()*np.sqrt(2))
    
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#    ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
    
    fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
    fig.savefig('flux_t='+str(t)+' density='+str(rho)+' E='+str(E)+'.eps') # Only to be saved as either .png or .jpg


def multiple(wdir):
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
        Xmesh, Ymesh = np.meshgrid(D.x1.T,D.x2.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.vx1.T,newdims,method='linear')
        velycong = T.congrid(D.vx2.T,newdims,method='linear')
        plt.gca().quiver(xcong, ycong, velxcong, velycong,color='w')
        plt.savefig('20_'+str(i)+'.png') # Only to be saved as either .png or .jpg
        plt.close()
#    show()
def temp(wdir):
    D = pp.pload(30,w_dir=wdir)
    
    I = pp.Image()
    flux=D.prs*1.67e-7/D.rho/1000000/1.3806488e-23
#    print flux.shape
#    flux= (flux-np.mean(flux))*5+np.mean(flux)*5.3
#    flux=nd.gaussian_filter(flux,sigma=(4,4),order=0)
    I.pldisplay(D, flux,x1=D.x1, \
                x2=D.x2,label1='l offset (pc)',label2='b offset (pc)',                                    \
                title='Temperature',                        
                cbar=(True,'vertical'))
#    savefig('MHD_Blast.png') # Only to be saved as either .png or .jpg
    plt.show()

def td(ty,t,E,rho,sigma,wdir):
    D = pp.pload(t,w_dir=wdir)   
    if ty == 'flux':
        for k in range(D.rho.shape[2]):
            for j in range(D.rho.shape[1]):
                for i in range(D.rho.shape[0]):
    #                if (i-D.rho.shape[0]/2)**2+(j-D.rho.shape[1]/2)**2+(k-D.rho.shape[2]/2)**2 > 30**2:
    #                    D.rho[k,j,i] = rho
                    if D.rho[k,j,i] > 10:
                        D.rho[k,j,i] = rho
        
        flux = D.rho*(D.bx1**2+D.bx2**2+D.bx3**2)**1.25*(D.vx1**2+D.vx2**2+D.vx3**2)**0
        flux = (flux-np.mean(flux))*1+np.mean(flux)*1
        flux = flux.sum(axis=1)
        flux = nd.gaussian_filter(flux,sigma=(sigma,sigma),order=0)
        
        fig = plt.figure(figsize=(6.3,5.5))
        ax = AxesGrid(fig, 111,  # similar to subplot(122)
                nrows_ncols=(1, 1),
                axes_pad=0.0,
                label_mode="1",
                share_all=True,
                cbar_location="right",
                cbar_mode="edge",
                cbar_size="2%",
                cbar_pad="0%",
                )
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
        neg = ax[0].imshow(flux.T*10000,origin='lower',extent=[D.x2[0],D.x2[-1],D.x3[0],D.x3[-1]])
        levels = np.arange(-5.5, -4, 0.5)
#        plt.contour(np.log10(flux).T, levels,
#                     origin='lower',
#                     linewidths=2,
#                     extent=[D.x1[0],D.x1[-1],D.x2[0],D.x2[-1]])
#        cmaps = [plt.get_cmap("spring"), plt.get_cmap("winter")]
#    
#        im = grid[0].imshow(Z, extent=extent, interpolation="nearest",
#                            cmap=cmaps[0])
    
        ax.cbar_axes[0].colorbar(neg)
    
        for cax in ax.cbar_axes:
            cax.toggle_label(True)
            cax.axis[cax.orientation].set_label(r'log(S)')
            
#        cbar = fig.colorbar(neg,ax=ax)
#        cbar.set_label(r'log(S)')
        ax.axes_llc.set_xlabel('l offset (pc)')
        ax.axes_llc.set_ylabel('b offset (pc)')
#        ax.axes_llc.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.99,bottom=0.08,left=0.11,right=0.91)
        fig.savefig('t'+str(t)+'_density'+str(rho)+'_E'+str(E)+'.eps')
        plt.show()
#        fig.clear()
#        ax = fig.add_subplot(111)
#        ax.plot((np.rot90(np.eye(len(flux)))*flux.T).sum(axis=0))
#        ax.set_xlim(0,256)
#        fig.savefig('change.eps')
        
    elif ty == 'rho':
#        t = 850
        fig = plt.figure(figsize=(7,6))
        ax = fig.add_subplot(111)
        neg = ax.imshow(np.log10(D.rho[:,:,64]).T,origin='lower',extent=[D.x2[0],D.x2[-1],D.x3[0],D.x3[-1]])
        cbar = fig.colorbar(neg,ax=ax)
        cbar.set_label(r'log($\mathregular{\rho/cm^{-3}}$)')
        ax.set_xlabel('x offset (pc)')
        ax.set_ylabel('y offset (pc)')
        ax.set_title(r't='+str(t)+r'$\ \mathregular{\rho}$='+str(rho)+' E='+str(E))
        fig.subplots_adjust(top=0.9,bottom=0.1,left=0.11,right=0.97)
        T = pp.Tools()
        newdims = 2*(20,)
        Xmesh, Ymesh = np.meshgrid(D.x2.T,D.x3.T)
        xcong = T.congrid(Xmesh,newdims,method='linear')
        ycong = T.congrid(Ymesh,newdims,method='linear')
        velxcong = T.congrid(D.bx1[:,:,64].T,newdims,method='linear')
        velycong = T.congrid(D.bx2[:,:,64].T,newdims,method='linear')
        plt.gca().quiver(xcong, ycong, velxcong, velycong,color='w')
        plt.show()
        fig.savefig('rho-t='+str(t)+'_density='+str(rho)+'_E='+str(E)+'.eps') # Only to be saved as either .png or .jpg
    #    close()
    else: 
        print(D.x1.shape)
#        arr = np.meshgrid(D.x1,D.x2,D.x3)
#        mlab.points3d(arr[0][0:256:8,0:256:8,0:256:8], arr[1][0:256:8,0:256:8,0:256:8], arr[2][0:256:8,0:256:8,0:256:8], D.rho[0:256:8,0:256:8,0:256:8])
        vol = mlab.pipeline.volume(mlab.pipeline.scalar_field((D.bx1**2+D.bx2**2+D.bx3**2)*D.rho))
        # ctf = ColorTransferFunction()
        # ctf.add_hsv_point(-8, 0.8, 1, 1)
        # ctf.add_hsv_point(-6.5, 0.45, 1, 1)
        # ctf.add_hsv_point(-5.4, 0.15, 1, 1)

        # vol._volume_property.set_color(ctf)
        # vol._ctf = ctf
        # vol.update_ctf = True
        # otf = PiecewiseFunction()

        # otf.add_point(-8, 0)
        # otf.add_point(-5.7, 0.082)
        # otf.add_point(-5.4, 0.0)

        # vol._otf = otf
        # vol._volume_property.set_scalar_opacity(otf)
#        mlab.contour3d(D.prs)
    #    mlab.quiver3d(D.bx1, D.bx2, D.bx3)
    #    src = mlab.pipeline.vector_field(D.bx1, D.bx2, D.bx3)
    #    mlab.pipeline.vectors(src, mask_points=20000, scale_factor=30.)
    #    mlab.outline()
#        mlab.savefig(str(t)+'.obj')
    
#==================================main========================================  
if __name__=='__main__':   
    #font = {'family' : 'serif',  
    #        'weight' : 'normal',  
    #        'size'   : 12,  
    #        }  
    
    choose = 'single' #single or multiple or temp
    ty     = 'flux'    
    t      = 5 
    E      = 1.3
    rho    = 1.0
    sigma  = 1
    
    wdir = './'
    nlinf = pp.nlast_info(w_dir=wdir)
#    single(ty,t,E,rho,sigma,wdir)
#    analyze(ty,t,E,rho,sigma,wdir)
    td(ty,t,E,rho,sigma,wdir)
