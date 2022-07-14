# gas_lift.py - A collection of functions to simulate gas lift control approaches
# GMX, 14 July 2022
#
# This file contains a collection of update and output functions for gas lifted
# oil wells to be used with the control package.

def glow_update(t, x, u, params={}):
  """Gas-lifted oil well dynamics
    Parameters
    ----------
    x: array
         System state: x1, mass of gas in the annulus [kg]
                       x2, mass of gas in the tubing [kg]
                       x3, mass of liquid in the tubing [kg]
        
    u: array
         System input: u1, opening of gas-lift choke [-]
                       u2, opening of production choke [-]
                       Pres, reservoir pressure [bar]
                       Pgs, gas lift source pressure [bar]
  
    Returns
    -------
    y: array
         The time derivatives of x1, x2 and x3
    """
  
  import numpy as np
  
  # Set up the system parameters

  R       = params.get('R',      8314.0)      # 
  g       = params.get('g',         9.81)     # 
  mu      = params.get('mu',        3.64e-03) # 
  rhoL    = params.get('rhoL',    760.0)      # 
  MG      = params.get('MG',       16.7)      # 
  Ta      = params.get('Ta',      348.0)      # 
  Va      = params.get ('Va',      64.34)     # 
  La      = params.get('La',     2048.0)      # 
  Pgs     = params.get('Pgs',     140.0)
  Vt      = params.get('Vt',       25.03)
  Sbh     = params.get('Sbh',       0.0314)
  Lbh     = params.get('Lbh',      75.0)
  Tt      = params.get('Tt',      369.4)
  GOR     = params.get('GOR',     0.0)
  Pres    = params.get('Pres',  160.0)
  wbarres = params.get('wbarres',18.0)
  Dt      = params.get('Dt',      0.134)
  Lt      = params.get('Lt',   2048.0)
  PI      = params.get('PI',      0.134)
  Kgs     = params.get('Kgs',     9.98e-05) 
  Kinj    = params.get('Kgs',     1.40e-04)
  Kpr     = params.get('Kgs',     2.90e-03)
  epsilon = params.get('epsilon', 1.5e-06)

  x1 = x[0]   # mass of gas in the annulus [kg]
  x2 = x[1]   # mass of gas in the tubing [kg]
  x3 = x[2]   # mass of liquid in the tubing [kg]

  u1   = u[0] # opening of gas-lift choke [-]
  u2   = u[1] # opening of production choke [-]
  Pres = u[2] # reservoir pressure [bar]
  Pgs  = u[3] # gas lift source pressure [bar]
  P0   = u[4] # pressure after the choke valve [bar]

  # Define auxiliary equations
  Pat = R*Ta*x1/MG/Va                            #(4)
  Pab = Pat + x1*g*La/Va                         #(5)
  rhoGab = Pab*MG/R/Ta                           #(6)
  rhoGin = Pgs*MG/R/Ta                           #(7)
  wGin = Kgs*u2*np.sqrt(rhoGin*max(Pgs - Pat,0)) #(8)
  rhoGt = x2/(Vt + Sbh*Lbh - x3/rhoL)            #(9)
  Ptt = rhoGt*R*Tt/MG                            #(10)
  rhobarmix = (x2 + x3 - rhoL*Sbh*Lbh)/Vt        #(11)
  alphabarL = (x3 - rhoL*Sbh*Lbh)/(Vt*rhoL)      #(12)
  alphamGb = GOR/(GOR + 1)                               #(13)
  Ubarslt = 4*(1 - alphamGb)*wbarres/rhoL/np.pi/Dt**2    #(14)
  Ubarsgt = 4(wGin + alphamGb*wbarres)/rhoGt/np.pi/Dt**2 #(15)
  Ubarmt = Ubarslt + Ubarsgt                             #(16)
  Ret = rhobarmix*Ubarmt*Dt/mu                           #(17)
  lambdat = fluids.friction.Haaland(Ret, epsilon/Dt)     #(18)
  Ft = alphabarL*lambdat*rhobarmix*Ubarmt**2*Lt/2/Dt     #(19)
  Ptb = Ptt + rhobarmix*g*Lt + Ft                        #(20)
  wGinj = Kinj*np.sqrt(rhoGab*max(Pab - Ptb,0))          #(21)
  Ubarlb = wbarres/(rhoL*Sbh)                            #(22)
  Reb = rhoL*Ubarlb*Db/mu                                #(23) 
  lambdab = fluids.friction.Haaland(Reb, epsilon/Db)     #(24)
  Fb = lambdab*rhoL*Ubarlb**2*Lbh/2/Db                   #(25)
  Pbh = Ptb + Fb + rhoL*g*Lbh                            #(26)
  wres = PI*max(Pres - Pbh,0)                            #(27)
  wLres = (1 - alphamGb)*wres                            #(28)
  wGres = alphamGb*wres                                  #(29)
  rhoGtb = Ptb*MG/R/Tt                                   #(30)
  alphaLb = wLres*rhoGtb/(wLres*rhoGtb + (wGinj + wGres)*rhoL)        #(31)
  alphaLt = 2*alphabarL - alphaLb                                     #(32)
  rhomixt = alphaLt + (1 - alphaLt)*rhoGt                             #(33)
  wout = Kpr*u1*np.sqrt(rhomixt*max(Ptt - P0,0))                      #(34)
  Qout = wout/rhomixt                                                 #(35)
  alphamGt = (1 - alphaLt)*rhoGt/(alphaLt*rhoL + (1 - alphaLt)*rhoGt) #(36)
  wGout = alphamGt*wout                                               #(37)
  wLout = (1 - alphamGt)*wout                                         #(38)
  
  # Define the ODEs
  dx1dt = wGin - WGinj                           #(1)
  dx2dt = WGinj + WGres - WGout                  #(2)
  dx3dt = wLres - WLout                          #(3)
 
  return [dx1dt, dx2dt, dx3dt]
