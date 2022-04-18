import matplotlib.pyplot as plt
from control.matlab import linspace, tf, feedback, lsim
from pandas import DataFrame
from numpy import interp
    
def pidtest(Kp, Ki, Kd, 
            N     = 10,
            num   = [1e-1],
            den   = [20**2, 2*0.5*20, 1],
            theta = 0,
            tpar  = [0, 10*20, 1000, 25], #start, stop, N, t0
            test  = True):
    r'''Simulates a servo loop with a plant and a PID controller (filtered derivative) and an open 
    loop for process reaction curve data.

    Parameters
    ----------
    Kp : float
         Proportional gain, [-] 
    Ki : float
         Integral gain, [1/time]   
    Kd : float
         Derivative gain, [time] 
    N  : float
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
        Gc  = tf([Kd, Kp, Ki], [Kd/Kp/N, 1, 0])
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
