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
    """
    Loads the initial conditions at the core of the star 

    Inputs:
    m (float): Enclosed mass at the core 
    comp (Tuple; (float, float float)): composition of the star in (X, Y, Z) 
    Pc (float): Initial guess for pressure at the core of the star 
    Tc (float): Initial guess for temperature at the core of the star 

    returns: 
    Tuple (float, float, float, float): Initial conditions at the core of the star in the format (rc, Pr, Lr, Tr)
    
    """
    # assigns variables
    X, Y, Z = comp   
    gamma2 = 5/3
    nabla_ad= 0.4

    mu = 4/(3+5*X) # calculate mean molecular weight from composotion    
    rhoc = rho_eq(Pc, Tc, comp) #calculate central density 

    kappac = calc_k(np.log10(Tc), np.log10(rhoc)) # determine opacity 
    
    rc = (3/(4*np.pi*rhoc))**(1/3) * m**(1/3) # calcualte radius 

    # calc energy gen rate for luminosity 
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
        Tr = T_rad # if radiative, set temp as radiative temperature 
        nabla = nabla_rad
    else:
        Tr = T_conv # "", convective 
        nabla = nabla_ad 
        
    return (rc, Pr, Lr, Tr)


def load2(M, comp, Ls, Rs):
    """
    Loads the initial conditions at the surface of the star. 

    Inputs:
    M (float): Enclosed mass at the surface (total star mass) 
    comp (Tuple; (float, float float)): composition of the star in (X, Y, Z) 
    Ls (float): Initial guess for luminosity at the surface of the star 
    Rs (float): Initial guess for the total radius of the star 

    returns: 
    Tuple (float, float, float, float): Initial conditions at the surface of the star in the format (Rs, Ps, Ls, Ts)
    
    """
    
    X, Y, Z = comp
    mu = 4/(3+5*X-Z)
    Ts = (Ls/(4*np.pi*(Rs**2)*sb))**(1/4) 

    # ensure that temperature is within range of opacity table 
    if not (3.75 <= np.log10(Ts) <= 7.5): 
        print(Ts)
        print("breaks at lum")  
        return -np.inf

    # optimize for surface pressure, density 
    def P_opt(rhos, Ts):
        """
        Optimizes the pressure as a function of density and temperature at the surface of the star using a system of two equations for surface pressure. 

        Inputs: 
        rhos (float): density at the surface of the star; independent variable of iteration in fsolve
        Ts (float): Temperature at the surface of the star 

        returns (float): 1-(P2/P1); factor to minimize within fsolve to get P2 and P1 to match, zero for P2/P1 = 1
        """

        # calculate rho, opacity table variables 
        rhos = rhos[0]
        R = rhos / (Ts*1e-6)**3 # calc R from T and rho 
        logR = np.log10(R) # take log of R 
        logT=np.log10(Ts)

        # checks conditions for opacity to make sure optimizer does not go in the wrong direction 
        if not(-8.0 < logR < 1.0):
            print("breaks at logR in P_opt")
            return np.inf
        if not (3.75 <= logT <= 7.5): 
            print("breaks at logT in P_opt")
            return np.inf
        # calc opacity, pressure with two equations 
        kappa = calc_k(np.log10(Ts), np.log10(rhos))
        P1 = (G*M)/(Rs**2) * (2/3) * (1/kappa) # pressure at surface of star 
        P2 = (rhos*Na*k*Ts)/mu + (1/3)*a*Ts**4
        return 1-(P1/P2)

    # optimize pressure as a function of rho at the surface of the star 
    rho_opt = fsolve(P_opt, x0 = [1.0*1e-10], args=(Ts))[0]
    kappa = calc_k(np.log10(Ts), np.log10(rho_opt))

    # calc pressure using one of the equations 
    Ps = (rho_opt*Na*k*Ts/mu) + (a*Ts**4)/3
    
    return (Rs, Ps, Ls, Ts)
    
    