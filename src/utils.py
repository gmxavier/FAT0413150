import matplotlib.pyplot as plt
from control.matlab import linspace, tf, feedback, lsim
from pandas import DataFrame
from numpy import interp
    
def pidtest(Kp, Ki, Kd, 
            alpha = 0.02,
            num   = [1e-1],
            den   = [20**2, 2*0.5*20, 1],
            theta = 0,
            tpar  = [0, 10*20, 1000, 25], #start, stop, N, t0
            test  = True):
    r'''Simulates a servo loop with a plant and an ideal PID controller and an open 
    loop for process reaction curve data.

    Parameters
    ----------
    Kp : float
         Proportional gain, [-] 
    Ki : float
         Integral gain, [1/time]   
    Kd : float
         Derivative gain, [time] 
    alpha: float
         Derivative filter factor [-]
    num : array
         Plant transfer function numerator coefficients        
    den : array
         Plant transfer function denominator coefficients        
    theta : float
         Plant dead time, [time]
    tpar : array
         Time parameters (start time, stop time, samples number and step time)
    test : boolean
         If test = True, it returns the servo control response.
         If test = False, it returns the process reaction curve.

    Returns
    -------
    retval : data frame
         Process reaction curve data (Time, Input and Output)

    Notes
    -----
  

    Example
    --------

    >>> pidtest(100,5,2000)


    Reference
    ----------

    '''            
    Gp = Gcl = tf(num, den)
    plt.title('Process reaction curve')
    if test == True:
        Gc  = tf([Kd, Kp, Ki], [alpha*Kd/Kp, 1, 0])
        Gcl = feedback(Gc*Gp, 1)
        plt.title('Servo control response \n Kp = {0} | Ki = {1} | Kd = {2}'.format(Kp,Ki,Kd))
    time = linspace(tpar[0], tpar[1], tpar[2])
    u = ud = 0*time
    u[time>tpar[3]] = 1
    ud[time>(tpar[3]+theta)] = 1
    y = lsim(Gcl, ud, time)[0]
    plt.plot(time, y, time, u)
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.legend(['Output', 'Input'])
    if test != True:
        retval = DataFrame({'Time': time, 'Input': u, 'Output': y})
        return retval

def callender(K, tau, theta, 
              type_of_plant='FODT',
              type_of_control='regulatory', 
              type_of_controller='PI', 
              criteria=1):
    r'''Returns the PI controller parameters from the rule of Callender et al. (1935/6).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)
    criteria : integer
         If criteria = 1, the rule gives decay ratio = 0.015 and period of decaying 
         oscillation = 5.10*theta.
         If criteria = 2, the rule gives decay ratio = 0.043 and period of decaying 
         oscillation = 6.28*theta.

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]

    Notes
    -----
    Applicable to theta/tau = 0.3.

    Example
    --------

    >>> callender(K=1, tau=10, theta=3)
    [0.18933333333333333, 0.01733821733821734]

    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''
    if criteria not in [1,2]:
        return []
    if criteria == 1:
        Kp = 0.568/(K*theta)
        Ki = Kp/(3.64*theta)
        return [Kp, Ki]
    if criteria == 2:
        if type_of_controller == 'PI':
            Kp = 0.690/(K*theta)
            Ki = Kp/(2.45*theta)
            return [Kp, Ki]
        if type_of_controller == 'PID':
            Kp = 1.066/(K*theta)
            Ki = Kp/(1.418*theta)
            Kd = Kp*(0.47*theta)
            return [Kp, Ki, Kd]

def ziegler_nichols(K, tau, theta, 
                   type_of_plant='FODT',
                   type_of_control='regulatory', 
                   type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rule of Ziegler and Nichols (1942).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----
    Applicable to theta/tau <= 1.0.

    Example
    --------

    >>> Kp, Ki, Kd = ziegler_nichols(K=1.25, tau=4, theta=0.9, type_of_controller='PID'); [Kp, Kp/Ki, Kd/Kp]
    [4.266666666666667, 1.8, 0.45]
    
    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''
    if type_of_controller == 'P':
        Kp = (1/K)*(tau/theta)
        return [Kp]
    if type_of_controller == 'PI':
        Kp = 0.9*(1/K)*(tau/theta)
        Ki = Kp/(3.3*theta)
        return [Kp, Ki]
    if type_of_controller == 'PID':
        Kp = 1.2*(1/K)*(tau/theta)
        Ki = Kp/(2*theta)
        Kd = Kp*0.5*theta
        return [Kp, Ki, Kd]    
 
def hazebroek_vanderwaerden(K, tau, theta, 
                            type_of_plant='FODT',
                            type_of_control='regulatory', 
                            type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rule of Hazebroek and Van der Waerden (1950).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----


    Example
    --------

    >>> hazebroek_vanderwaerden(K=1, tau=10, theta=3)
    [2.3333333333333335, 0.16339869281045755]

    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''    
    thetaovertau = theta/tau
    if thetaovertau < 0.2:
        return [] 
    if (thetaovertau >= 0.2) & (thetaovertau <= 3.4):    
        thetaovertau_ = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,
                         1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4]
        x1_ = [0.68,0.70,0.72,0.74,0.76,0.79,0.81,0.84,0.87,0.90,0.93,0.96,0.99,
               1.02,1.06,1.09,1.13,1.17,1.20,1.28,1.36,1.45,1.53,1.62,1.71,1.81]
        x2_ = [7.14,4.76,3.70,3.03,2.50,2.17,1.92,1.75,1.61,1.49,1.41,1.32,1.25,
               1.19,1.14,1.10,1.06,1.03,1.00,0.95,0.91,0.88,0.85,0.83,0.81,0.80]
        x1 = interp(thetaovertau, thetaovertau_, x1_)
        x2 = interp(thetaovertau, thetaovertau_, x2_)
        Kp = x1/(K*thetaovertau)
        Ki = Kp/(x2*theta)
        return [Kp, Ki]        
    if thetaovertau > 3.4:
        Kp = 1/(K*thetaovertau)*(0.5*thetaovertau + 1)
        Ki = Kp/(theta/(1.6*theta - 1.2*tau))
        return [Kp, Ki]  
 
def oppelt(K, tau, theta, 
           type_of_plant='FODT',
           type_of_control='regulatory', 
           type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rule of Oppelt (1951).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----
    This function calculates the recommended value for Kp.

    Example
    --------

    >>> oppelt(K=1, tau=10, theta=3)
    [1.5666666666666669, 0.15729585006693445]

    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''    
    thetaovertau = theta/tau
    if thetaovertau < 1:
        Kp = (1/K)*(0.77/thetaovertau - 1)
        Ki = Kp/(3.32*theta)
        return [Kp, Ki]  
    if thetaovertau > 1:
        Kp = (1/K)*(0.77/thetaovertau - 1)
        Ki = Kp/(1.66*theta)
        return [Kp, Ki]  

def moros_oppelt(K, tau, theta, 
                 type_of_plant='FODT',
                 type_of_control='regulatory', 
                 type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rules of Moros (1999).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----


    Example
    --------

    >>> moros_oppelt(K=1, tau=10, theta=3)


    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''    
    thetaovertau = theta/tau
    Kp = 0.8/(K*thetaovertau)
    Ki = Kp/(3*theta)
    return [Kp, Ki]  

def moros_rosenberg(K, tau, theta, 
                    type_of_plant='FODT',
                    type_of_control='regulatory', 
                    type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rules of Moros (1999).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----


    Example
    --------

    >>> moros_rosenberg(K=1, tau=10, theta=3)


    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''    
    thetaovertau = theta/tau
    Kp = 0.91/(K*thetaovertau)
    Ki = Kp/(3.3*theta)
    return [Kp, Ki] 

def cohen_coon(K, tau, theta, 
               type_of_plant='FODT',
               type_of_control='regulatory', 
               type_of_controller='PI'):
    r'''Returns the PI controller parameters from the rule of Cohen and Coon (1953).

    Parameters
    ----------
    K : float
         Static gain of the process reaction curve, [-] 
    tau : float
         Time constant (lag) of the process reaction curve, [time]
    theta : float
         Dead time of the process reaction curve, [time]
    type_of_plant : string
         Type of the plant model   
    type_of_control : string
         Type of the control loop (regulatory or servo)
    type_of_controller : string
         Type of the controller (P, PI, PD, PID)

    Returns
    -------
    Kp : float
         Proportional gain, [-]

    Ki : float
         Integral gain, [1/time]
         
    Kd : float
         Derivative gain, [time]         

    Notes
    -----
    Applicable to theta/tau <= 1.0.

    Example
    --------

    >>> Kp, Ki, Kd = cohen_coon(K=1.25, tau=4, theta=0.9, type_of_controller='PID'); [Kp, Kp/Ki, Kd/Kp]
    [4.266666666666667, 1.8, 0.45]
    
    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''
    if type_of_controller == 'P':
        Kp = (1/K)*(tau/theta + 1/3)
        return [Kp]
    if type_of_controller == 'PI':
        Kp = (1/K)*(0.9*tau/theta + 1/12)
        Ki = Kp/((theta*(30 + 3*(theta/tau))/(theta*(9 + 20*(theta/tau)))))
        return [Kp, Ki]
    if type_of_controller == 'PID':
        Kp = (1/K)*((4/3)*(tau/theta) + 1/4)
        Ki = Kp/((theta*(32 + 6*(theta/tau))/(theta*(13 + 8*(theta/tau)))))
        Kd = Kp*(4/(theta*(11 + 2*(theta/tau))))
        return [Kp, Ki, Kd]    
