import numpy as np
from rho import * 
from energies import * 
from constants import *
from star_params import * 
from nabla_check import check_nabla

def dPdMr_eq(Mr, r):
    """
    Calculates the pressure derivative dPdMr with respect to mass

    Inputs:
    Mr (float): Mass 
    r (float): radius 

    returns: (float) dPdMr
    """
    return (-G*Mr)/(4*np.pi*(r**4))
    
def drdMr_eq(Mr, r, P, T): 
    """
    Calculates the radius derivative drdMr with respect to mass

    Inputs:
    Mr (float): Mass 
    r (float): radius 
    P (float): Pressure
    T (float): Temperature

    returns: (float) drdMr
    """
    
    rho = rho_eq(P, T, comp) 
    return 1/(4*np.pi*(r**2) * rho)
    
def dLrdMr_eq(Mr, P, T):
    """
    Calculates the luminosity derivative dLrdMr with respect to mass

    Inputs:
    Mr (float): Mass 
    P (float): Pressure
    T (float): Temperature

    returns: (float) dLrdMr
    """
    rho = rho_eq(P, T, comp)
    X, Y, Z = comp
    e_CNO_calc = e_CNO(X, Y, Z, rho, T)
    e_PP_calc = e_PP(X, Y, Z, rho, T)
    
    return e_CNO_calc + e_PP_calc
    
def dTdMr_eq(Mr, T, r, P, L):
    """
    Calculates the temperature derivative dTdMr with respect to mass

    Inputs:
    Mr (float): Mass 
    T (float): Temperature
    r (float): radius 
    P (float): Pressure
    L (float) Luminosity

    returns: (float) dLdMr
    """
    nabla = check_nabla(comp, Mr, T, P, L) # uses check_nabla function to get stellar gradient if radiative or convective  
    return -(G*Mr*T*nabla)/(4*np.pi*(r**4)*P)

def derivs(m, params):
    """
    Calculates and compiles the four coupled differential equations for dPdMr, drdMr, dLrdMr, and dTdMr at an enclosed mass m having corresponding params r, P, L T
    
    Inputs: 
    m (float): enclosed mass coordinate
    params (list): in the form [r, P, L, T], corresponding to the radius, pressure, luminosity, and temperature at the current mass coordiante m 

    returns: (tuple) of all calculated derivatives in the form (drdMr, dPdMr, dLrdMr, dTdMr)
    """
    
    r, P, L, T = params 
    dPdMr = dPdMr_eq(m, r)
    drdMr = drdMr_eq(m, r, P, T)
    dLrdMr = dLrdMr_eq(m, P, T)
    dTdMr = dTdMr_eq(m, T, r, P, L)
    
    return (drdMr, dPdMr, dLrdMr, dTdMr)