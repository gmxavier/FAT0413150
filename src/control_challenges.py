
# control_challenges.py - A collection of functions to simulate examples of
# some levels of challenges at https://janismac.github.io/ControlChallenges/
# GMX, 1 September 2024
#
# This file contains a collection of update and output functions to be used 
# with the control package.

try: 
  import control as ct
except:
  import pip
  package_names=['control'] #packages to install
  pip.main(['install'] + package_names)
  import control as ct


import matplotlib.pyplot as plt        
  
def tutorial_update(t, x, u, params):
  import numpy as np
  """Tutorials dynamics
  Parameters
  ----------
  x: array
       System state: p_x     position [m]
                     v_x     speed    [m/s]
        
  u: array
       System input: F_t     force    [N]
  
  Returns
  -------
  y: array
       The time derivatives of p_x, v_x
  """   
  
  # Set up the system parameters
  g       = params.get('g',         9.81)     # gravitational constant   [m/s2]
  mu      = params.get('mu',         0.0)     # coefficient of friction   [-]
  m_t     = params.get('m_t',        1.0)     # block mass   [kg]
  theta   = params.get('theta',      0.0)     # plane angle   [rad]
  
  p_x     = x[0] # position [m]
  v_x     = x[1] # speed    [m/s]

  F_t = u[0]     # force    [N]
 
  # Define the ODEs
  p_xdot  = v_x
  v_xdot  = (F_t - m_t*g*np.sin(theta) - mu*m_t*g*np.cos(theta))/m_t
 
  return np.array([p_xdot, v_xdot])

def tutorial_control_output(t, x, u, params):
  import numpy as np
  """Tutorial generic PD controller
  Parameters
  ---------        
  u: array
      System input: SP   setpoint                      [-]
                    PV   process variable              [-]
                    dPV  time der. of process variable [-]
  y: array
      System output: CO  controller output             [-]
  
  Returns
  -------
  y: array
      The value of CO.
    """

  # Set up the system parameters
  K_P = params.get('K_P', 1.0e+00) # proportional gain [-]
  K_D = params.get('K_D', 1.0e+00) # derivative gain   [s]
 
  SP  = u[0] # setpoint                      [-]
  PV  = u[1] # process variable              [-]
  dPV = u[2] # time der. of process variable [-]

  # Define the auxiliary equations
  E  = SP - PV # error              [-]
  dE = -dPV    # time der, of error [-]

  # Define the outputs
  CO = K_P*E + K_D*dE # controller output [-] 
 
  return np.array([CO])

def tutorial_plot(sys, T, u, x0, params):
  import numpy as np
  try: 
    import control as ct
  except:
    import pip
    package_names=['control'] #packages to install
    pip.main(['install'] + package_names)
  import control as ct
  
  t, y = ct.input_output_response(sys, T, u, x0, params)
  
  xlabel = 'Time $t$ [s]'
  ylabel = ['Force $F_t$ [N]',
            'Position $p_x$ [m]',
            'Speed $v_x$ [m/s]'
            ]

  res = {xlabel : t}

  fig, axes = plt.subplots(3, 1, figsize=(22, 22))
  axe = axes.ravel()

  for k, ax in enumerate(axe):
    scale = 1
    ax.plot(t, y[k]/scale)
    res[ylabel[k]] = y[k]/scale
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel[k])
    ax.grid(True, linestyle='dotted')

  return res
