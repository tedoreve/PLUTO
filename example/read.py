import numpy as np
#b=np.fromfile("rho0.dbl",dtype=np.float)
#np.savetxt('xx.txt',b) 
a=np.random.uniform(0,1,size=16384)
a.tofile("rho0.dbl")

#---------------------------------------------------------------------
b=np.linspace(-0.5,0.5,129)
c=np.linspace(1,128,len(b)-1)
d=np.array([c,b]).T
f=open('grid0.out','w')
f.write('# GEOMETRY:   CARTESIAN\n')
f.write(str(len(b)-1)+'\n')
for i in range(len(c)):
    f.write(str(int(c[i]))+'  '+str(b[i])+'  '+str(b[i+1])+'\n')
f.write(str(len(b)-1)+'\n')
for i in range(len(c)):
    f.write(str(int(c[i]))+'  '+str(b[i])+'  '+str(b[i+1])+'\n')
f.write('1\n')
f.write('1 0.0 1.0')
f.close()
#print d
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#Spectral_Index=-5/3
#Minimum_Lambda = 4
#Norm      = 0.1
#IR_Cutoff = 2.0
#NXLim     = 2*(int(len(c)/2/Minimum_Lambda))
#NYLim     = 2*(int(len(c)/2/Minimum_Lambda))
#NXLim_2   = NXLim/2
#NYLim_2   = NYLim/2
#
#Amplitude = np.zeros([NXLim + 1, NYLim + 1])
#Phase     = np.zeros([NXLim + 1, NYLim + 1])
#twoPI     = np.pi*2
#x_end     = 1
#y_end     = 1
#for i in range(NXLim):
#    for j in range(NYLim):
#        x1r  = np.random.uniform(0,1,size=1)
#        x2r  = np.random.uniform(0,1,size=1)
#        xdbl = np.sqrt(-2*np.log(x1r))*np.cos(twoPI*x2r)
#        Phase[i,j]=twoPI*np.random.uniform(0,1,size=1)
#        scrh1= IR_Cutoff**2+(twoPI*(i-NXLim_2)/x_end)**2+(twoPI*(i-NYLim_2)/y_end)**2
#        Amplitude[i,j]=np.sqrt(Norm*scrh**(Spectral_Index/2))*xdbl
#for j2 in range():












