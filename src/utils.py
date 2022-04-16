import matplotlib.pyplot as plt
from control.matlab import linspace, tf, feedback, lsim
from pandas import DataFrame
from google.colab import data_table
data_table.enable_dataframe_formatter()

def pidtest(Kp, Ki, Kd, 
            alpha = 0.02,
            num   = [1e-1],
            den   = [20**2, 2*0.5*20, 1],
            theta = 0,
            t0 = 25,
            tpar  = [0, 10*20, 0.2],
            test = True):
  
  Gp = Gcl = tf(num, den)
  plt.title('Process reaction curve')
  if test == True:
    Gc  = tf([Kd, Kp, Ki], [alpha*Kd/Kp, 1, 0])
    Gcl = feedback(Gc*Gp, 1)
    plt.title('Servo control response \n Kp = {0} | Ki = {1} | Kd = {2}'.format(Kp,Ki,Kd))
  time = linspace(tpar[0], tpar[1], tpar[2])
  u = 0*time
  u[time>t0] = 1
  ud[time>(t0+theta)] = 1
  y = lsim(Gcl, ud, time)[0]
  plt.plot(time, y, time, u)
  plt.xlabel('Time')
  plt.ylabel('Amplitude')
  plt.legend(['Output', 'Input'])
  if test != True:
    d = DataFrame({'Time': time, 'Input': u, 'Output': y})
    return data_table.DataTable(d, include_index = False)
