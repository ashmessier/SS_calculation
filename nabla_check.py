from constants import * 
from rho import rho_eq
from opacities import calc_k
import numpy as np
from energies import e_CNO, e_PP


def check_nabla(comp, m, T, P, L, return_conv_rad = False):

    """
    Checks if the star is radiaitve or convective at a given mass coordiante m with corresponding temperature T, pressure P, and luminoisty L. Calculates the temperature
    corresponding to convective and radiative energy transport, as well as the gradients nabla_rad and nabla_ad. If nabla_ad is larger than nabla_rad, the star is radiative.
    If nabla_ad is larger than nabla_rad, the star is convective. 

    Inputs: 
    comp (tuple): composition of the star in (X, Y, Z)
    m (float): enclosed mass coordinate m 
    T (float): Temperature 
    P (float): Pressure
    L (float): Luminosity 
    return_conv_rad (Optional, bool) = False

    returns: 
    nabla (float): stellar gradient, either corresponding to nabla_rad or nabla_ad depending on the relative values of the two.
    if return_conv_rad = True; returns "radiative" if the star is radative at the mass m and "convective" if the star is convective
    
    """
    
    X, Y, Z = comp
    nabla_ad = 0.4 # from 1/(1-gamma); gamma = 5/3

    rho = rho_eq(P, T, comp) #calculate central density
    
    kappac = calc_k(np.log10(T), np.log10(rho)) # determine opacity 
    e_CNO_calc = e_CNO(X, Y, Z, rho, T) # CNO energy gen rate 
    e_PP_calc = e_PP(X, Y, Z, rho, T) # PP energy gen rate 
    
    ec = e_CNO_calc+e_PP_calc # determine total energy gen rate 

    # calculates radiative temperature, temp gradient 
    T_rad = (T**4 - (1/(2*a*c)) * (3/(4*np.pi))**(2/3) * kappac*ec*rho**(4/3)*m**(2/3))**(1/4)
    nabla_rad = 3/(16*np.pi*a*c) * (P*kappac)/(T_rad**4) * (L)/(G*m) 

    # calculates convective temperature 
    log_T_conv = np.log10(T) - (np.pi/6)**(1/3) * G * (nabla_ad * rho**(4/3))/(P) * (m**(2/3))
    T_conv = 10**(log_T_conv)

    Tr = 0 # initialize temperature output ; to be assigned to T_conv or T_rad later 
    nabla = 0 # initializes nabla output; ""
    
    if nabla_ad > nabla_rad: # checks for radiative condition
        Tr = T_rad # assigns Tr output 
        nabla = nabla_rad # assigns stellar gradient output 
        if return_conv_rad: # if return_conv_rad = True, returns string corresponding to condition 
            return "radiative"
    else: # else, convective, same as above
        Tr = T_conv
        nabla = nabla_ad 
        if return_conv_rad: 
                return "convective"

    return nabla