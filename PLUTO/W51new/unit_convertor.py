from astropy import units as un
from astropy import constants as con
import numpy as np
import matplotlib.pyplot as plt
#          定义的参量               数值备注       真正在程序中用的值       CGS单位               简略单位
##define  UNIT_DENSITY            1 CONST_mp           1.67e-24           g     cm^-3     
##define  UNIT_LENGTH             1 CONST_pc           3.09e18                  cm
##define  UNIT_VELOCITY          10000 km/s            1.00e9                   cm      s^-1
#         UNIT_B                                       4.58e-3            g^0.5 cm^-0.5 s^-1      G
#         UNIT_t                 98 years              3.09e9                           s
#         UNIT_P                                       1.67e-6            g     cm^-1   s^-2      Ba=0.1Pa
#         UNIT_e                                       4.93e49            g     cm^2    s^-2      erg
#         UNIT_M                                       4.93e31            g
UNIT_DENSITY = 1*un.M_p
UNIT_LENGTH  = 1*un.pc
UNIT_VELOCITY= 1e4*un.km/un.s
UNIT_B = 4.58e3*un.uG
UNIT_t = 98*un.yr
UNIT_P = 1.67e-6*un.Ba
UNIT_E = 4.93e49*un.erg
UNIT_M = 2.48e-2*con.M_sun

n   = 1*un.M_p
l   = 1*un.pc
v   = 6000*un.km/un.s
B   = 9*un.uG
t   = 100*un.yr
P   = 1*un.Ba
E   = 1e51*un.erg
M   = 14*con.M_sun

n   /= UNIT_DENSITY
l   /= UNIT_LENGTH
v   /= UNIT_VELOCITY
B   /= UNIT_B
t   /= UNIT_t
P   /= UNIT_P
E   /= UNIT_E
M   /= UNIT_M

#print('n = ', n, '\n'
#      'l = ', l, '\n'
#      'v = ', v, '\n'
#      'B = ', B, '\n'
#      't = ', t, '\n'
#      'P = ', P, '\n'
#      'E = ', E, '\n'
#      'M = ', M, '\n')

def g(x,u,o1,o2):
    return np.exp(-0.5*((x-u)/((np.arctan(x-u)+np.pi)/2/np.pi*(o2-o1)+o1))**2)
x=np.linspace(1,300,100)
o1=50
o2=1
u=200
plt.plot(x,g(x,u,o1,o2))  