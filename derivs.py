import numpy as np
from rho import * 
from energies import * 
from constants import *
from star_params import * #params0, Ms, M_intersect
from nabla_check import check_nabla

def dPdMr_eq(Mr, r):
    return (-G*Mr)/(4*np.pi*(r**4))
    
def drdMr_eq(Mr, r, P, T): # make sure this is ok, something is turning negative 
    rho = rho_eq(P, T, comp) 
  #  print("rho in drdMR eq,", rho)
    return 1/(4*np.pi*(r**2) * rho)
    
def dLrdMr_eq(Mr, P, T): # get this guy from somewhere (I think last hw)
    rho = rho_eq(P, T, comp)
    X, Y, Z = comp
    e_CNO_calc = e_CNO(X, Y, Z, rho, T)
    e_PP_calc = e_PP(X, Y, Z, rho, T)
  #  print(e_CNO_calc+e_PP_calc)
    
    return e_CNO_calc + e_PP_calc
    
def dTdMr_eq(Mr, T, r, P, L): # huh 
    nabla = check_nabla(comp, Mr, T, P, L)
   # nabla=0.4
    return -(G*Mr*T*nabla)/(4*np.pi*(r**4)*P)

def derivs(m, params):
    
    r, P, L, T = params # make sure calls everything 
    dPdMr = dPdMr_eq(m, r)
    drdMr = drdMr_eq(m, r, P, T)
    dLrdMr = dLrdMr_eq(m, P, T)
    dTdMr = dTdMr_eq(m, T, r, P, L)

    #print("derivs", drdMr, dPdMr, dLrdMr, dTdMr)
    
    return (drdMr, dPdMr, dLrdMr, dTdMr)