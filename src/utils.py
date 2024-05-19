import matplotlib.pyplot as plt
from control.matlab import linspace, tf, feedback, lsim
from pandas import DataFrame
from numpy import interp
    
def pidtest(Kp, Ki, Kd, 
            a     = 0.1,
            M     = 1,
            num   = [1e-1],
            den   = [20**2, 2*0.5*20, 1],
            theta = 0,
            tpar  = [0, 10*20, 1000, 25], #start, stop, N, t0
            test  = True,
            out   = False):
    r'''Simulates a closed-loop (servo control) with a PID controller (filtered derivative) or an open-loop
    step test for process reaction curve data.

    Parameters
    ----------
    Kp : float
         Proportional gain, [-] 
    Ki : float
         Integral gain, [1/time]   
    Kd : float
         Derivative gain, [time] 
    a  : float
         Derivative filter factor [-]
    M  : float
         Step amplitude
    num : array
         Process reaction curve transfer function numerator coefficients        
    den : array
         Process reaction curve transfer function denominator coefficients        
    theta : float
         Process reaction curve transfer function dead time, [time]
    tpar : array
         Time parameters (start time, stop time, samples number and step time)
    test : boolean
         If test = True, it returns the closed-loop (servo control) response.
         If test = False, it returns the open-loop step test response.
    out  : boolean
         If out = True, it returns the time, the input and the output.

    Returns
    -------
    retval : data frame
         Closed-loop (servo control) or process reaction curve response 
         (Time, Input and Output)

    Notes
    -----
  

    Example
    --------

    >>> pidtest(100,5,2000)


    Reference
    ----------

    '''            
    Gp = Gcl = tf(num, den)
    plt.title('Open-loop step test response')
    if test == True:
        Gc  = tf([Kd, Kp, Ki], [a*Kd/Kp, 1, 0])
        Gcl = feedback(Gc*Gp, 1)
        plt.title('Closed-loop (servo control) response \n Kp = {0} | Ki = {1} | Kd = {2}'.format(Kp,Ki,Kd))
    time = linspace(tpar[0], tpar[1], tpar[2])
    u = 0*time
    u[time>tpar[3]] = M
    ud = 0*time
    ud[time>(tpar[3]+theta)] = M
    y = lsim(Gcl, ud, time)[0]
    plt.plot(time, y, time, u)
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.legend(['Controller output', 'Sensor output'])
    if test == True:
        plt.legend(['Sensor output', 'Internal setpoint'])
    if out == True:
        retval = DataFrame({'Time': time, 'Input': u, 'Output': y})
        return retval

def fom(inpval):
    r''' Calculates several figures of merit of a closed-loop control (e.g. overshoot, decay ratio, 
    integral of absolute error, integral of squared error, integral of time-weighted absolute error, 
    integral of time-weighted squared error).

    Parameters
    ----------
    inpval : data frame
         Closed-loop (servo control) response (Input is the ISP) 
         (Time, Input and Output)

    Returns
    -------
    retval : tuple
         Figures of merit of a closed-loop (servo control)
         (OS, DR, IAE, ISE, ITAE, ITSE)

    Notes
    -----
  

    Example
    --------

    >>> inpval = pidtest(100,5,2000, out=True)
    >>> fom(inval)


    Reference
    ----------

    '''            
    t = inpval['Time']
    u = inpval['Input']
    y = inpval['Output']
    threshold = u.iloc[-1]
    if u.iloc[-1] > 0:
        aux = np.where(np.where([(y - np.roll(y,1) > 0) & (y - np.roll(y,-1) > 0)],y, 0)> threshold, y,np.nan)
        OS = aux[0] - threshold
        DR = (aux[1] - threshold)/OS
    if u.iloc[-1] < 0:
        aux = np.where(np.where([(y - np.roll(y,1) < 0) & (y - np.roll(y,-1) < 0)],y, 0)< threshold, y,np.nan)
        OS = aux[0] - threshold
        DR = (aux[1] - threshold)/OS    
    #OS = min(y - u.iloc[-1]) if u.iloc[-1]<0 else max(y - u.iloc[-1])
    #DR = min(y[y<OS] - u.iloc[-1]) if u.iloc[-1]<0 else max(y[y<OS] - u.iloc[-1])
    IAE = sum(abs(u - y))*(max(t)-min(t))/len(t)
    ISE = sum((u - y)**2)*(max(t)-min(t))/len(t)
    ITAE = sum(abs(u - y)*t)*(max(t)-min(t))/len(t)
    ITSE = sum((u - y)**2*t)*(max(t)-min(t))/len(t)
    retval = (OS,DR,IAE,ISE,ITAE,ITSE)
    return retval
