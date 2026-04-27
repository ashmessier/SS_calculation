from opacities import *
from energies import * 
from rho import * 
from constants import *


import numpy as np
import scipy 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp



def load1(m, comp, Pc, Tc):
    
    X, Y, Z = comp
    # Initial guesses
   
    gamma2 = 5/3
    nabla_ad= (gamma2 - 1)/gamma2

    mu = 4/(3+5*X) # calculate mean molecular weight from composotion    
    rhoc = rho_eq(Pc, Tc, comp) #calculate central density 

    kappac = calc_k(np.log10(Tc), np.log10(rhoc)) # determine opacity 
    
    rc = (3/(4*np.pi*rhoc))**(1/3) * m**(1/3) # this was in my notes from last day before BD lecture 
    
    e_CNO_calc = e_CNO(X, Y, Z, rhoc, Tc)
    e_PP_calc = e_PP(X, Y, Z, rhoc, Tc)
                
    ec = e_CNO_calc + e_PP_calc

    Lr = ec * m # calc luminosity from energy gen rate 
    Pr = Pc - (3*G)/(8*np.pi) * (4*np.pi/3 * rhoc)**(4/3) * m**(2/3) # calc pressure 
    
    # determine if convective or adiabatic depending on nabla_ad vs nabla_rad 
    T_rad = (Tc**4 - (1/(2*a*c)) * (3/(4*np.pi))**(2/3) * kappac*ec*rhoc**(4/3)*m**(2/3))**(1/4)
    nabla_rad = 3/(16*np.pi*a*c) * (Pr*kappac)/(T_rad**4) * (Lr)/(G*m) 

    log_T_conv = np.log10(Tc) - (np.pi/6)**(1/3) * G * (nabla_ad * rhoc**(4/3))/(Pc) * (m**(2/3))
    T_conv = 10**(log_T_conv)

    Tr = 0 # initialize temperature 
    nabla = 0 
    
    if nabla_ad > nabla_rad: 
        Tr = T_rad # if 
        nabla = nabla_rad
    else:
        Tr = T_conv
        nabla = nabla_ad 

  #  print("Rc", rc)
        
    return (rc, Pr, Lr, Tr)


def load2(M, comp, Ls, Rs):
    X, Y, Z = comp
    mu = 4/(3+5*X-Z)
   # print("in load2 Ls:", Ls, "Rs",Rs)
    Ts = (Ls/(4*np.pi*(Rs**2)*sb))**(1/4) # from luminosity guess?
  #  print(np.log10(Ts))
   # print("TS in load2", Ts)
    if not (3.75 <= np.log10(Ts) <= 7.5): 
        print(Ts)
        print("breaks at lum")
        return -np.inf

    def P_opt(rhos, Ts):
        rhos = rhos[0]
        R = rhos / (Ts*1e-6)**3 # calc R from T and rho 
        logR = np.log10(R) # take log of R 
        logT=np.log10(Ts)
     #   print("logR", logR)

        if not(-8.0 < logR < 1.0):
            print("breaks at logR in P_opt")
            return np.inf
        if not (3.75 <= logT <= 7.5): 
            print("breaks at logT in P_opt")
            return np.inf
        
        kappa = calc_k(np.log10(Ts), np.log10(rhos))
        P1 = (G*M)/(Rs**2) * (2/3) * (1/kappa) # pressure at surface of star 
        P2 = (rhos*Na*k*Ts)/mu + (1/3)*a*Ts**4
        return 1-(P1/P2)

    rho_opt = fsolve(P_opt, x0 = [1.0*1e-10], args=(Ts))[0]
    kappa = calc_k(np.log10(Ts), np.log10(rho_opt))
   # Ps = (G*M)/(Rs**2) * (2/3) * (1/kappa) # pressure at surface of star  # why is this negative 
    Ps = (rho_opt*Na*k*Ts/mu) + (a*Ts**4)/3
  #  print("pressure in load2", Ps)
   # print("load2 output",Rs, Ps, Ls, Ts)
    
    return (Rs, Ps, Ls, Ts)
    
    