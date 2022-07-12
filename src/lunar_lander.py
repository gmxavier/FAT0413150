try: 
  import control as ct
except:
  import pip
  package_names=['control'] #packages to install
  pip.main(['install'] + package_names)
  import control as ct


import matplotlib.pyplot as plt        
  
def lunar_update(t, x, u, params={}):
  import numpy as np
  """Lunar lander dynamics
  Parameters
  ----------
  x: array
       System state: p_x     horizontal position [m]
                     v_x     horizontal speed    [m/s]
                     p_y     vertical position   [m]
                     v_y     vertical speed      [m/s]
                     theta   heading             [rad]
                     v_theta heading speed       [rad/s]
                     m_t     total mass          [kg]
                     m_f     fuel mass           [kg]
        
  u: array
       System input: F_t   main engine thrust    [-]
                     F_l   lateral engine thrust [-]
  
  Returns
  -------
  y: array
       The time derivatives of p_x, v_x, p_y, v_y, theta, v_theta, m_t and m_f
  """   
  
  # Set up the system parameters
  # Source: https://web.aeromech.usyd.edu.au/AMME3500/Course_documents/material/tutorials/Assignment%204%20Lunar%20Lander%20Solution.pdf
  g       = params.get('g',        1.6)     # moon gravitational constant   [m/s2]
  L       = params.get('L',        4.0)     # distance to center of gravity [m]
  J       = params.get('J',      100.0e+03) # moment of inertia             [kg m2]
  I_sp    = params.get('I_sp',     3.0e+03) # thruster specific impulse     [N/(kg/s)]
  F_t_max = params.get('F_t_max', 44.0e+03) # main engine max thrust        [N]
  F_l_max = params.get('F_l_max',  0.5e+03) # lateral engine max thrust     [N]

  #p_x_0   =   0.5e+03 # initial horizontal position  [m]
  #v_x_0   =   0.0e+03 # initial horizontal speed     [m/s]
  #p_y_0   = 160.0e+03 # initial vertical position    [m]
  #v_y_0   =  -0.7e+03 # initial vertical speed       [m/s]
  #theta   =   0.0e+03 # initial heading              [rad]
  #v_theta =   0.0e+03 # initial heading speed        [rad/s]
  #m_t_0   =  15.0e+03 # initial total mass           [kg]
  #m_f_0   =   8.0e+03 # initial fuel mass            [kg]

  p_x     = x[0] # horizontal position [m]
  v_x     = x[1] # horizontal speed    [m/s]
  p_y     = x[2] # vertical position   [m]
  v_y     = x[3] # vertical speed      [m/s]
  theta   = x[4] # heading             [rad]
  v_theta = x[5] # heading speed       [rad/s]
  m_t     = x[6] # total mass          [kg]
  m_f     = x[7] # fuel mass           [kg]
 
  F_t = u[0]*(np.sign(m_f) + 1)/2 # main engine thrust    [-]
  F_l = u[1]*(np.sign(m_f) + 1)/2 # lateral engine thrust [-]

  # Define the auxiliary equations
  c_f = ((F_t + np.abs(F_l))/I_sp)

  # Define the ODEs
  p_xdot     = v_x
  v_xdot     = (F_l*np.cos(theta) - F_t*np.sin(theta))/m_t
  p_ydot     = v_y
  v_ydot     = (F_l*np.sin(theta) + F_t*np.cos(theta))/m_t - g
  thetadot   = v_theta
  v_thetadot = (L/J)*F_l
  m_tdot     = -c_f
  m_fdot     = -c_f
 
  return np.array([p_xdot, v_xdot, 
                   p_ydot*(np.sign(p_y) + 1)/2, # set vertical speed to zero when on ground
                   v_ydot*(np.sign(p_y) + 1)/2, # set vertical acceleration to zero when on ground
                   thetadot, v_thetadot, m_tdot, m_fdot])

def lunar_engine_output(t, x, u, params={}):
  import numpy as np
  """Lunar lander engine
  Parameters
  ---------        
  u: array
      System input: T_t   main engine throtlle    [-]
                    T_l   lateral engine throtlle [-]
                    m_f   fuel mass               [kg]

  y: array
      System output: F_t   main engine thrust    [N]
                     F_l   lateral engine thrust [N]                       
  
  Returns
  -------
  y: array
      The values of F_t and F_l
    """

  # Set up the system parameters
  # Source: https://web.aeromech.usyd.edu.au/AMME3500/Course_documents/material/tutorials/Assignment%204%20Lunar%20Lander%20Solution.pdf
  F_t_max = params.get('F_t_max', 44.0e+03) # main engine max thrust        [N]
  F_l_max = params.get('F_l_max',  0.5e+03) # lateral engine max thrust     [N]
 
  T_t = u[0] # main engine throtlle    [-]
  T_l = u[1] # lateral engine throtlle [-]
  m_f = u[2] # fuel mass               [kg]

  # Define the auxiliary equations
  T_t = np.clip(T_t,  0, 1)*(np.sign(m_f) + 1)/2
  T_l = np.clip(T_l, -1, 1)*(np.sign(m_f) + 1)/2 
  
  # Define the outputs
  F_t = T_t*F_t_max # main engine thrust        [N]
  F_l = T_l*F_l_max # lateral engine thrust     [N]
 
  return np.array([F_t, F_l])

def lunar_control_output(t, x, u, params={}):
  import numpy as np
  """Lunar lander generic PD controller
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

def landing_plot(sys, T, u, x0, params={}):
  import numpy as np
  try: 
    import control as ct
  except:
    import pip
    package_names=['control'] #packages to install
    pip.main(['install'] + package_names)
  import control as ct
  
  v_y_max = params.get('v_y_max', 2.5) # max vertical speed [m/s]  
  
  t, y = ct.input_output_response(sys, T, u, x0)
  
  xlabel = 'Time $t$ [s]'
  ylabel = ['Main engine throttle $T_t$ [%]',
            'Lateral engine throttle $T_l$ [%]',
            'Horizontal position $p_x$ [km]',
            'Horizontal speed $v_x$ [m/s]',
            'Vertical position $p_y$ [km]',
            'Vertical speed $v_y$ [m/s]',
            'Heading $\\theta$ [rad]',
            'Heading speed $v_\\theta$ [rad/s]',
            'Total mass $m_t$ [Mg]',
            'Fuel mass $m_f$ [Mg]'
            ]

  res = {xlabel : t}

  fig, axes = plt.subplots(5, 2, figsize=(22, 22))
  axe = axes.ravel()

  for k, ax in enumerate(axe):
    scale = 1e3 if k in [2,4,8,9] else 1e-2 if k in [0,1] else 1
    y[k] = np.clip(y[k], 0, 1) if k in [0,1] else y[k]
    if ((k == 0) & (y[4][-1] < 0) & (-y[5][-1] < v_y_max)):
      ax.set_title('Congratulations, the lunar lander touched down at {:.1f} m/s!'.format(-y[5][-1]))
    elif ((k == 0) & (y[4][-1] < 0) & (-y[5][-1] > v_y_max)):
      ax.set_title('Sorry, the lunar lander crashed at {:.1f} m/s!'.format(-y[5][-1]))
    elif ((k == 0) & (y[4][-1] > 0)):
      ax.set_title('Sorry, the lunar lander did not touch down! Actually it is now at {:.3f} km of altitude!'.format(y[4][-1]/1e3))
    ax.plot(t, y[k]/scale)
    res[ylabel[k]] = y[k]/scale
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel[k])
    ax.grid(True, linestyle='dotted')

  return res
