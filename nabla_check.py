from constants import * 
from rho import rho_eq
from opacities import calc_k
import numpy as np
from energies import e_CNO, e_PP
#from star_params import comp


def check_nabla(comp, m, T, P, L, return_conv_rad = False):
   # print("nabla for T:", T)
    X, Y, Z = comp
    nabla_ad = 0.4

    rho = rho_eq(P, T, comp) #calculate central density
    
    kappac = calc_k(np.log10(T), np.log10(rho)) # determine opacity 
    e_CNO_calc = e_CNO(X, Y, Z, rho, T)
    e_PP_calc = e_PP(X, Y, Z, rho, T)
    
    ec = e_CNO_calc+e_PP_calc
    
    T_rad = (T**4 - (1/(2*a*c)) * (3/(4*np.pi))**(2/3) * kappac*ec*rho**(4/3)*m**(2/3))**(1/4)
    nabla_rad = 3/(16*np.pi*a*c) * (P*kappac)/(T_rad**4) * (L)/(G*m) 

    log_T_conv = np.log10(T) - (np.pi/6)**(1/3) * G * (nabla_ad * rho**(4/3))/(P) * (m**(2/3))
    T_conv = 10**(log_T_conv)

    Tr = 0 # initialize temperature 
    nabla = 0
    
    if nabla_ad > nabla_rad: 
        Tr = T_rad # if 
        nabla = nabla_rad
        if return_conv_rad: 
            return "radiative"
       # print("radiative")
    else:
        Tr = T_conv
        nabla = nabla_ad 
       # print("convective")
        if return_conv_rad: 
                return "convective"

    return nabla