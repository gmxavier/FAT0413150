# dynamic_models.py - A collection of dynamic models examples
# GMX, 16 August 2021
#
# This file contains a collection of full nonlinear models of several different
# processes.

import numpy as np
try:
  import fluids
except:
  !pip install fluids
  import fluids.friction

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

def cstr_update(t, x, u, params={}):
  """CSTR dynamics

    Parameters
    ----------
    x: array
         System state: CA, molar concentration of component A at the reactor
                       in kmol/m3
                       T, reactor temperature in K
        
    u: array
         System input: F, feed flowrate in m3/h
                       CA0, molar concentration of component A at the feed
                       in kmol/m3
                       T0, feed tempeature in K
                       Tc, reactor cooling jacket temperature in K
  

    Returns
    -------
    y:  array
        The time derivatives of CA and T

    """

  # Set up the system parameters

  Ea  = params.get('Ea',11850.0)   # activation energy [kcal/kmol] 
  R   = params.get('R',1.98589)    # ideal gas constant [kcal/kmol/K]
  k_0 = params.get('k_0',3.5e+07)  # Arrhenius constant[1/h]
  dH  = params.get('dH',-5960.0)   # reaction enthalpy [kcal/kmol]
  rho = params.get('rho',1000.0)   # density [kg/m3]
  Cp  = params.get('Cp',0.480)     # heat capacity [kcal/(kg*K)]
  UA  = params.get ('UA',145.0)    # global heat transfer coefficient times heat transfer area [kcal/(K*h)]
  V   = params.get('V',1.0)        # volume [m3]
 
  CA= x[0]                         # molar concentration of component A at the feed [kmol/m3]
  T = x[1]                         # reactor temperature [K]

  F   = u[0]                       # feed flowrate [m3/h]
  CA0 = u[1]                       # molar concentration of component A at the feed, [kmol/m3]
  T0  = u[2]                       # reactor temperature [K]
  Tc  = u[3]                       # reactor cooling jacket temperature, [K]


  # Define the ODEs
  dCAdt = ((F/V)*(CA0-CA)) - (k_0*(CA*np.exp(-Ea/(R*T))))
  dTdt = ((F/V)*(T0-T))-((UA)/(V*Cp*rho))*(T-Tc)+((-dH*V*k_0)/(Cp*rho)*((np.exp(-Ea/(R*T))*CA)))
 
  return [dCAdt, dTdt]

def berg_update(t, x, u, params={}):
  """Insulin dynamics based on Bergman minimal model.

    Parameters
    ----------
    x : array
         System state: G, the deviation of the blood glucose concentration from
                       basal levels measured in mg dL−1
                       I, the deviation of the blood insulin concentration from 
                       basal levels measured in mU L−1
                       X, a proportionality variable that describes the insulin 
                       concentration in a remote compartment measured in mU L−1
    u : array
         System input: U, the insulin input from an external source and into the
                       bloodstream, measured in mU min−1
                       D, a disturbance variable representing the intake of 
                       glucose from external sources, typically food, measured 
                       in mg dL−1min−1

    Returns
    -------
    y : array
        The time derivatives of G, I and X

    """

  # Set up the system parameters
  Gb = params.get('Gb', 81.)      # mg/dL
  Ib = params.get('Ib', 15.)      # mU/L
  P1 = params.get('P1', 0.0287)   # 1/min
  P1_s = params.get('P1_s', 0.02) # 1/min
  P2 = params.get('P2', 0.0287)   # 1/min 
  P3 = params.get('P3', 5.035e-5) # L/mU/min^2  
  n = params.get('n', 5/54)       # 1/min  
  V1 = params.get('V1', 12)       # L  

  # Define variables for patient state and inputs
  G = x[0]                           # blood glucose concentration
  I = x[1]                           # blood insulin concentration
  X = x[2]                           # insulin concentration in a remote 
                                     # compartment
  U = u[0]                           # insulin input
  D = u[1]                           # intake of glucose

  # Define the ODEs
  dGdt = -P1*G - X*(G+Gb) + D
  dIdt = -n*(I+Ib) + U/V1
  dXdt = -P2*X + P3*I

  return [dGdt, dIdt, dXdt]

def polyeth_update(t, x, u, params={}):
  """Gas phase polyethylene reactor dynamics.

    Dadebo, S. A., Bell, M. L., McLellan, P. J., & McAuley, K. B. (1997). 
    Temperature control of industrial gas phase polyethylene reactors. 
    Journal of Process Control, 7(2), 83–95.
    https://www.sciencedirect.com/science/article/pii/S0959152496000169

    McAuley, K. B. (1991). Modelling, Estimation and Control of Product 
    Properties in a Gas Phase Polyethylene Reactor. Ph.D. Thesis.
    https://macsphere.mcmaster.ca/handle/11375/8939

    Parameters
    ----------
    x : array
         System state: In_con, molar concentration of inert components in the 
                       gas phase in mol/L
                       M1_con, molar concentration of ethylene in the 
                       gas phase in mol/L
                       Y1, moles of active site type 1, mol
                       Y2, moles of active site type 2, mol
                       T, reactor temperature in K
                       Tw1, temperature of cooling water stream leaving stage 1
                       of heat exchanger in K
                       Tg1, temperature of recycle gas stream leaving stage 1 
                       of heat exchanger in K
    u : array
         System input: Fc, flow rate of catalyst in kg/s
                       Tfeed, feed temperature in K
                       Twi, cooling water temperature in K                       

    Returns
    -------
    y : array
        The time derivatives of In_con, M1_con, Y1, Y2, T, Tw1 and Tg1

    """
  
  # Set up the system parameters
  Vg    = params.get('Vg', 500.)                   # m3
  Vp    = params.get('Vp',   0.5)                  # 
  Pv    = params.get('Pv',  17.)                   # atm
  Bw    = params.get('Bw', 7.0e4)                  # kg
  kp0   = params.get('kp0', 85.0e-3)               # m3/(mol s) 
  Ea    = params.get('Ea', 9e4*4.1868)             # J/mol  
  Cpm1  = params.get('Cpm1', 11*4.1868)            # J/(mol K)  
  Cv    = params.get('Cv', 7.5)                    # atm^-0.5
  Cpw   = params.get('Cpw', 4.1868e3)              # J/(kg K)
  CpIn  = params.get('CpIn', 6.9*4.1868)           # J/(kg K)
  Cppol = params.get('Cppol', 0.85e3*4.1868)       # J/(kg K)
  kd1   = params.get('kd1', 1e-4)                  # s^-1
  kd2   = params.get('kd2', 1e-4)                  # s^-1
  Mw1   = params.get('Mw1', 28.05e-3)              # kg/mol
  Mw    = params.get('Mw', 3.314e4)                # kg
  Mg    = params.get('Mg', 6060.5)                 # mol
  MrCpr = params.get('MrCpr', 1.4*4.1868e7)        # J/K
  Hreac = params.get('Hreac', -894*4.1868e3)       # J/kg
  UA    = params.get('UA', 1.14*4.1868e6)          # J/(K s)
  FIn   = params.get('FIn', 5.)                    # mol/s
  FM1   = params.get('FM1', 190.)                  # mol/s
  Fg    = params.get('Fg', 8500.)                  # mol/s
  Fw    = params.get('Fw', 3.11e5*18e-3)           # kg/s
  Tf    = params.get('Tf', 360.)                   # K
  RR    = params.get('RR', 8.20575e-5)             # (m3 atm)/(mol K)
  R     = params.get('R', 8.314)                   # J/(mol K)
  ac    = params.get('ac', 0.548)                  # mol/kg

  # Define variables for reactor state and inputs
  In_con = x[0]     # molar concentration of inert components in the gas phase
  M1_con = x[1]     # molar concentration of ethylene in the gas phase
  Y1     = x[2]     # moles of active site type 1
  Y2     = x[3]     # moles of active site type 2
  T      = x[4]     # reactor temperature
  Tw1    = x[5]     # temperature of cooling water stream leaving stage 1 of heat exchanger
  Tg1    = x[6]     # temperature of recycle gas stream leaving stage 1 of heat exchanger

  Fc     = u[0]     # flow rate of catalyst in kg/s
  Tfeed  = u[1]     # feed temperature in K
  Twi    = u[3]     # cooling water temperature in K

  # Define the algebraic equations
  bt   = Vp * Cv * np.sqrt((M1_con+In_con) * RR * T - Pv)
  RM1  = M1_con * kp0 * np. exp(-Ea/R*(1/T-1/Tf)) * (Y1+Y2)
  Cpg  = M1_con/(M1_con + In_con) * Cpm1 + In_con/(M1_con + In_con) * CpIn
  Hf   = FM1 * Cpm1 * ( Tfeed - Tf) + FIn * CpIn * (Tfeed - Tf)
  Hg1  = Fg * (Tg1 - Tf) * Cpg
  Hg0  = (Fg + bt) * (T - Tf) * Cpg
  Hr   = Hreac * Mw1 * RM1
  Hpol = Cppol * (T - Tf) * RM1 * Mw1

  # Define the ODEs
  dIn_condt = (FIn - In_con/(M1_con + In_con) * bt)/Vg
  dM1_condt = (FM1 - M1_con/(M1_con + In_con) * bt - RM1)/Vg
  dY1dt     = Fc * ac - kd1 * Y1 - RM1 * Mw1 * Y1/ Bw
  dY2dt     = Fc * ac - kd2 * Y2 - RM1 * Mw1 * Y2/ Bw
  dTdt      = (Hf + Hg1 - Hg0 - Hr - Hpol)/(MrCpr + Bw * Cppol)
  dTw1dt    = Fw/Mw * (Twi - Tw1) - UA/(Mw * Cpw) * (Tw1 - Tg1)
  dTg1dt    = Fg/Mg * (T - Tg1)   + UA/(Mg * Cpg) * (Tw1 - Tg1)

  return [dIn_condt, dM1_condt, dY1dt, dY2dt, dTdt, dTw1dt, dTg1dt]

def fourtank_update(t, x, u, params={}):
  """Four tank system dynamics.

    Ashok Kumar, B., Jeyabharathi, R., Surendhar, S., Senthilrani, S., & 
    Gayathri, S. (2019). Control of Four Tank System Using Model Predictive 
    Controller. 2019 IEEE International Conference on System, Computation, 
    Automation and Networking (ICSCAN).
    https://ieeexplore.ieee.org/abstract/document/8878700

    Parameters
    ----------
    x : array
         System state: h1, height of tank 1 in cm
                       h2, height of tank 2 in cm
                       h3, height of tank 3 in cm
                       h4, height of tank 4 in cm
    u : array
         System input: v1, controller signal to pump 1 in V
                       v2, controller signal to pump 2 in V

    Returns
    -------
    y : array
        The time derivatives of h1, h2, h3 and h4

    """
  
  # Set up the system parameters
  A1     = params.get('A1', 28.)                   # cm2
  A2     = params.get('A2', 32.)                   # cm2
  A3     = params.get('A3', 28.)                   # cm2
  A4     = params.get('A4', 32.)                   # cm2
  a1     = params.get('a1', 0.071)                 # cm2
  a2     = params.get('a2', 0.057)                 # cm2
  a3     = params.get('a3', 0.057)                 # cm2
  a4     = params.get('a4', 0.071)                 # cm2
  k1     = params.get('k1', 3.33)                  # (cm3/s)/V
  k2     = params.get('k2', 3.35)                  # (cm3/s)/V
  gamma1 = params.get('gamma1', 0.70)              # -
  gamma2 = params.get('gamma2', 0.60)              # -
  g      = params.get('g', 981)                    # cm/s2

  # Define variables for tank state and inputs
  h1 = x[0]     # height of tank 1 in cm
  h2 = x[1]     # height of tank 2 in cm
  h3 = x[2]     # height of tank 3 in cm
  h4 = x[3]     # height of tank 4 in cm

  v1 = u[0]     # controller signal to pump 1 in V
  v2 = u[1]     # controller signal to pump 1 in V

  # Define the algebraic equations

  # Define the ODEs
  dh1dt = -(a1/A1)*np.sqrt(2*g*h1) + (a3/A1)*np.sqrt(2*g*h3) + (gamma1*k1/A1)*v1
  dh2dt = -(a2/A2)*np.sqrt(2*g*h2) + (a4/A2)*np.sqrt(2*g*h4) + (gamma2*k2/A2)*v2
  dh3dt = -(a3/A3)*np.sqrt(2*g*h3) + ((1 - gamma2)*k2/A3)*v2
  dh4dt = -(a4/A4)*np.sqrt(2*g*h4) + ((1 - gamma1)*k1/A4)*v1

  return [dh1dt, dh2dt, dh3dt, dh4dt]

def nuclear_update(t, x, u, params={}):
  """Nuclear reactor dynamics.

    Gábor, A., Fazekas, C., Szederkényi, G., & Hangos, K. M. (2011). Modeling and 
    Identification of a Nuclear Reactor with Temperature Effects and Xenon Poisoning. 
    European Journal of Control, 17(1), 104–115.
    https://doi.org/10.3166/ejc.17.104-115
    
    Parameters
    ----------
    x : array
         System state: N,  neutron concentration in %
                       C,  concentration of the delayed neutron emitting nuclei in % 
                       Tf, temperature of the fuel in °C
                       Tm, average temperature of the moderator in °C
                       X,  xenon concentration in 1/cm3
                       I,  iodine concentration in 1/cm3
                       
    u : array
         System input: z,   rod position in m
                       Tin, temperature of the water entering the reactor in °C

    Returns
    -------
    y : array
        The time derivatives of N, C, Tf, Tm, nX and nI

    """
  
  # Set up the system parameters
  phi0   = params.get('phi0', 1.3e13)       # initial equilibrium neutron flux in 1/(cm2 s)
  sigmaX = params.get('sigmaX', 2.805e-18)  # microscopic absorption cross section in cm2
  alphaf = params.get('alphaf', -5.362e-3)  # temperature coefficient of the fuel in $/°C
  alpham = params.get('alpham', -2.075e-2)  # temperature coefficient of the moderator in $/°C
  A1     = params.get('A1', 0.1056)         # parameter group in 1/s
  A3     = params.get('A3', 0.8757)         # parameter group in 1/s
  p0     = params.get('p0', 0.0401)         # rod parameter 0 in $
  p1     = params.get('p1', -0.44)          # rod parameter 1 in $/m
  p2     = params.get('p2', -0.966)         # rod parameter 2 in $/m2
  LAMBD  = params.get('LAMBD', 2.18e-5)     # average generation time in s
  SIGMAf = params.get('SIGMAf', 0.3358)     # macroscopic fission cross section in 1/cm
  lambdI = params.get('lambdI', 2.849e-5)   # decay constant of iodine in 1/s
  lambdX = params.get('lambdX', 2.150e-5)   # decay constant of xenon in 1/s
  lambdC = params.get('lambdC', 7.728e-2)   # decay constant of the delayed neutron 
                                            # emitting nuclei in 1/s
  beta   = params.get('beta', 0.0065)       # fraction of delayed neutron group
  YI     = params.get('YI', 6.38e-2)        # iodine yield 
  YX     = params.get('YI', 2.2e-3)         # xenon yield
  
  if (t==0):
    N0   = x[0] # neutron concentration at t=0 in %
    Tf0  = x[2] # temperature of the fuel at t=0 in °C
    Tm0  = x[3] # temperature of the water leaving the reactor at t=0 in °C
    X0   = x[4] # xenon concentration per macroscopic fission cross section at t=0 in 1/cm2
    Tin0 = u[1] # temperature of the water entering the reactor at t=0 in °C

  # Define variables for reactor state and inputs
  N  = x[0]     # neutron concentration in %
  C  = x[1]     # concentration of the delayed neutron emitting nuclei in % 
  Tf = x[2]     # temperature of the fuel in °C
  Tm = x[3]     # temperature of the water leaving the reactor in °C
  X  = x[4]     # xenon concentration per macroscopic fission cross section in 1/cm2
  I  = x[5]     # iodine concentration per macroscopic fission cross section in 1/cm2

  z   = u[0]     # rod position in m
  Tin = u[1]     # temperature of the water entering the reactor in °C

  # Define the algebraic equations
  rho = alphaf*(Tf - Tf0) + alpham*(Tm - Tm0) + p2*z**2 + p1*z + p0 - (sigmaX/beta/SIGMAf)*(X - X0)

  # Define the ODEs
  dNdt  = (beta/LAMBD)*(rho - 1)*N + (beta/LAMBD)*C
  dCdt  = lambdC*(N - C)
  dTfdt = -A1*(Tf - Tm) + A1*((Tf0 - Tm0)/N0)*N
  dTmdt = A3*(Tf - Tm) - A3*((Tf0 - Tm0)/(Tm0 - Tin0))*(Tm - Tin)
  dIdt  = YI*N/N0*phi0 - lambdI*I
  dXdt  = YX*N/N0*phi0 + lambdI*I - lambdX*X - sigmaX*N/N0*phi0
                             
  return [dNdt, dCdt, dTfdt, dTmdt, dIdt, dXdt]

# USEFUL REFERENCES
#
# Control structure design for stabilizing unstable gas-lift oil wells
# https://folk.ntnu.no/skoge/publications/2012/jahanshahi-adchem12-gaslift/0110.pdf
#
# Modeling and Identification of a Nuclear Reactor with Temperature Effects and Xenon Poisoning
# https://www.sciencedirect.com/science/article/abs/pii/S0947358011705726
#
# Improved temperature control of a PWR nuclear reactor using an LQG/LTR based controller
# https://ieeexplore.ieee.org/document/1178712
#
# A simple dynamic model of the primary circuit in VVER plants for controller design purposes
# https://www.sciencedirect.com/science/article/abs/pii/S0029549306006789?via%3Dihub
#
# Drum-boiler dynamics
# http://www.dei.unipd.it/~picci/Files/TEACHING/CONTROLLO_DEI_PROCESSI/DrumBoiler.pdf
