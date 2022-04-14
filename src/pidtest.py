import scipy as sp
import numpy as np
from control.matlab import *
import matplotlib.pyplot as plt
import pandas as pd
from google.colab import data_table
data_table.enable_dataframe_formatter()

def pidtest(Kp, Ki, Kd, 
            alpha = 0.02,
            num   = [1e-1],
            den   = [20**2, 2*0.5*20, 1],
            theta = 0,
            t0 = 25,
            T  = np.arange(0, 10*20, 0.2),
            test = True):
  
  Gp = Gcl = tf(num, den)
  plt.title('Process reaction curve')
  if test == True:
    Gc  = tf([Kd, Kp, Ki], [alpha*Kd/Kp, 1, 0])
    Gcl = feedback(Gc*Gp, 1)
    plt.title('Servo control response \n Kp = {0} | Ki = {1} | Kd = {2}'.format(Kp,Ki,Kd))
  u = np.ones(len(T))
  u[T<=t0] = 0
  ud = np.ones(len(T))
  ud[T<=(t0+theta)] = 0
  y = lsim(Gcl, ud, T)[0]
  plt.plot(T, y, T, u)
  plt.xlabel('Time')
  plt.ylabel('Amplitude')
  plt.legend(['Output', 'Input'])
  if test != True:
    d = pd.DataFrame({'Time': T, 'Input': u, 'Output': y})
    return data_table.DataTable(d, include_index = False)
