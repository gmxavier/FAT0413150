import matplotlib.pyplot as plt
from control.matlab import linspace, tf, feedback, lsim
from pandas import DataFrame
    
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
         If test = True, it simulates the servo loop response.
         If test = False, it simulates the process reaction curve.

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

def callender356(K, tau, theta, 
                 type_of_plant='FODT',
                 type_of_control='Regulatory', 
                 type_of_controller='PI', 
                 criteria=1):
    r'''Calculates the PI controller parameters for a FOLPD process reaction 
    curve using the rule of Callender et al. (1935/6).

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
         Type of the control loop (Regulatory or Servo)
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

    >>> callender356(K=1, tau=10, theta=3)
    [0.18933333333333333, 10.92]  

    Reference
    ----------
    .. [1] O’Dwyer, A. Handbook of PI and PID Controller Tuning Rules. London:
       Imperial College Press, 2009.
    '''
    if criteria not in [1,2]:
        return []
    # Rule 1 (Table 9, p. 30)
    if criteria == 1:
        Kp = 0.568/(K*theta)
        Ki = Kp/(3.64*theta)
        return [Kp, Ki]
    # Rule 2 (Table 9, p. 30)
    if criteria == 2:
        Kc = 0.690/(K*theta)
        Ki = Kp/(2.45*theta)
        return [Kp, Ki]
