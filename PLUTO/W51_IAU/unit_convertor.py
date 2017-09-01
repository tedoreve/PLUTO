from astropy import units as un
from astropy import constants as con
import numpy as np
#          定义的参量               数值备注       真正在程序中用的值       CGS单位               简略单位
##define  UNIT_DENSITY            1 CONST_mp           1.67e-24           g     cm^-3     
##define  UNIT_LENGTH             1 CONST_pc           3.09e18                  cm
##define  UNIT_VELOCITY          10000 km/s            1.00e9                   cm      s^-1
#         UNIT_B                                       4.58e-3            g^0.5 cm^-0.5 s^-1      G
#         UNIT_t                 98 years              3.09e9                           s
#         UNIT_P                                       1.67e-6            g     cm^-1   s^-2      Ba=0.1Pa
#         UNIT_e                                       4.93e49            g     cm^2    s^-2      erg
#         UNIT_M                                       4.93e31            g
UNIT_DENSITY = 1*con.m_p/un.cm**3
UNIT_LENGTH  = 1*un.pc
UNIT_VELOCITY= 1e4*un.km/un.s
UNIT_B = (UNIT_VELOCITY*np.sqrt(4*np.pi*UNIT_DENSITY)).to(un.g**0.5*un.cm**-0.5*un.s**-1).value*un.G
UNIT_t = UNIT_LENGTH/UNIT_VELOCITY
UNIT_P = UNIT_DENSITY*UNIT_VELOCITY**2
UNIT_M = UNIT_DENSITY*UNIT_LENGTH**3
UNIT_E = UNIT_M*UNIT_VELOCITY**2
#UNIT_B = ((UNIT_E/UNIT_LENGTH**3)**0.5).value*un.G
UNIT_NU= UNIT_P*UNIT_t

n   = 2*con.m_p/un.cm**3
l   = 4*un.pc
v   = 490*un.km/un.s
B   = 9*un.uG
t   = 1000*un.yr
P   = 1*un.Ba
E_th= 0.96*un.erg
E   = 2e51*un.erg
M   = 15*con.M_sun
nu  = 2*un.uPa*un.s

n   /= UNIT_DENSITY
l   /= UNIT_LENGTH
v   /= UNIT_VELOCITY
B   /= UNIT_B
t   /= UNIT_t
P   /= UNIT_P
E   /= UNIT_E
M   /= UNIT_M
nu  /= UNIT_NU

print('n = ', n.to('').value, '\n'
      'l = ', l.to('').value, '\n'
      'v = ', v.to('').value, '\n'
      'B = ', B.to('').value, '\n'
      't = ', t.to('').value, '\n'
      'P = ', P.to('').value, '\n'
      'E = ', E.to('').value, '\n'
      'M = ', M.to('').value, '\n'
      'nu = ', nu.to('').value, '\n'
      )

