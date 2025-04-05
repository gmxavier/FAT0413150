# beer_fermentor.py - A collection of functions to simulate beer
# fermentation process
# GMX, 04 April 2025
#
# This file contains a collection of update and output functions for a
# beer fermentor to be used with the control package.

try: 
  import control as ct
except:
  import pip
  package_names=['control'] #packages to install
  pip.main(['install'] + package_names)
  import control as ct


import matplotlib.pyplot as plt        
  
def beer_update(t, x, u, params):
  import numpy as np
  """Beer fermentor dynamics
  Parameters
  ----------
  x: array
       System state: X_A     active cells concentration [g/L]
                     X_L     latent cells concentration [g/L]
                     X_D     dead cells concentration [g/L]
                     S       sugar concentration [g/L]
                     EtOH    ethanol concentration [g/L]
                     DY      diacetyl concentration [ppm]
                     EA      ethyl acetate concentration [ppm]
        
  u: array
       System input: T       temperature [°C]
  
  Returns
  -------
  y: array
       The time derivatives of X_A, X_L, X_D, S, EtOH, DY, EA
  """   
  
  # Set up the system parameters
  # Source: https://web.aeromech.usyd.edu.au/AMME3500/Course_documents/material/tutorials/Assignment%204%20Lunar%20Lander%20Solution.pdf
  # g       = params.get('g',        1.6)     # moon gravitational constant   [m/s2]
  # L       = params.get('L',        4.0)     # distance to center of gravity [m]
  # J       = params.get('J',      100.0e+03) # moment of inertia             [kg m2]
  # I_sp    = params.get('I_sp',     3.0e+03) # thruster specific impulse     [N/(kg/s)]
  # F_t_max = params.get('F_t_max', 44.0e+03) # main engine max thrust        [N]
  # F_l_max = params.get('F_l_max',  0.5e+03) # lateral engine max thrust     [N]

  # X_A     =  # active cells concentration [g/L]
  # X_L     =  # latent cells concentration [g/L]
  # X_D     =  # dead cells concentration [g/L]
  # S       =  # sugar concentration [g/L]
  # EtOH    =  # ethanol concentration [g/L]
  # DY      =  # diacetyl concentration [ppm]
  # EA      =  # ethyl acetate concentration [ppm]

  X_A     = x[0] # active cells concentration [g/L]
  X_L     = x[1] # latent cells concentration [g/L]
  X_D     = x[2] # dead cells concentration [g/L]
  S       = x[3] # sugar concentration [g/L]
  EtOH    = x[4] # ethanol concentration [g/L]
  DY      = x[5] # diacetyl concentration [ppm]
  EA      = x[6] # ethyl acetate concentration [ppm]
 
  # F_t = u[0]*(np.sign(m_f) + 1)/2 # main engine thrust    [-]
  # F_l = u[1]*(np.sign(m_f) + 1)/2 # lateral engine thrust [-]

  # Define the auxiliary equations
  # Source: https://doi.org/10.1016/S1474-6670(17)44433-8
  mu_x_0 = np.exp(108.31 - 31934.09*T)
  mu_x = mu_x_0*S/(k_x + EtOH)


  # Define the ODEs
  X_Adot     = mu_x*X_A - mu_DT*X_A + mu_L*X_L
  X_Ldot     = -mu_L*X_L
  X_Ddot     = -mu_SD*X_D + mu_DT*X_A
  Sdot       = -mu_S*X_A
  EtOHdot    = f*mu_e*X_A
  DYdot      = mu_DY*S*X_A - mu_AB*DY*EtOH
  EAdot      = Y_EA*mu_x*XA
 
  return np.array([X_Adot, X_Ldot, X_Ddot, Sdot, thetadot, EtOHdot, DYdot, EAdot])