from load1load2 import * 
from scipy.integrate import solve_ivp
import numpy as np
from star_params import * #params0, Ms, M_intersect
from derivs import * 


def int_load1_outwards(Rs, Pc, Ls, Tc): 

    """
    Function to integrate outwards from the core using the four initial guesses for surface radius (Rs), core pressure (Pc), total luminosity (Ls), core temperature (Tc).
    Calculations only reference Pc and Tc. 

    Inputs: 
    Rs (float): Radius at the surface of the star (does not reference)
    Pc (float): Pressure at the core of the star 
    Ls (float): Luminosity at the surface of the star (does not reference) 
    Tc (float): Temperature at the core of the star 
    """
        
    params = load1(M_core_init, comp, Pc, Tc) # uses load1 to call initial core conditions 

    # integrats params outwards using the derivs function from the core to the intersection mass 0.2 Ms
    sol_outwards = solve_ivp(derivs, t_span = [M_core_init, M_intersect], y0 = params, method = "Radau", dense_output=True)
    return sol_outwards

def int_load2_inwards(Rs, Pc, Ls, Tc):
    """
    Function to integrate inwards from the surface using the four initial guesses for surface radius (Rs), core pressure (Pc), total luminosity (Ls), core temperature (Tc).
    Calculations only reference Ls and Rs. 

    Inputs: 
    Rs (float): Radius at the surface of the star 
    Pc (float): Pressure at the core of the star (does not reference)
    Ls (float): Luminosity at the surface of the star
    Tc (float): Temperature at the core of the star (does not reference) 
    """

    params = load2(Ms, comp, Ls, Rs) 
    
    # integrates params inwards with derivs function from surface to 0.2 Ms 
    sol_inwards = solve_ivp(derivs, t_span = [Ms, M_intersect], y0 = params, method = "Radau", dense_output=True)
    return sol_inwards

def difference(params):
    """
    Integrates inwards and outwards using int_load2_inwards and int_load1_outwards functions to the intersection mass; calculates the difference between the inwards and 
    outwards integrations at the intersection mass. Function is to be optimized over using fsolve. 

    Inputs: 
    params (list): list in the form of [Pc, Tc, Ls, Rs] 

    returns: 
    (list) in the form [difference_R, difference_P, difference_L, difference_T], corresponding to the normalized difference between the inner and outer integrations 
    
    """
        
    Pc, Tc, Ls, Rs = params

    results_in = int_load2_inwards(Rs, Pc, Ls, Tc).y

    results_out = int_load1_outwards(Rs, Pc, Ls, Tc).y

    if not (3.75 <= np.log10(results_in[-1][-1]) <= 7.5):
        print("log T bad in load1")
        return np.ones(4) * np.inf

    # check radius
    if results_in[0][-1] < 0: 
        print("bad radius in")
        return -np.ones(4) * np.inf
    if results_out[0][-1] < 0: 
        print("bad radius out")
        return -np.ones(4) * np.inf

    # check luminosity 
    if results_in[2][-1] < 0: 
        print(results_in[2][-1])
        print("bad lum in")
        return -np.ones(4) * np.inf
    if results_out[2][-1] < 0: 
        print("bad lum out")
        return -np.ones(4) * np.inf

            # check temp 
    if results_in[3][-1] < 0: 
        print(results_in[3][-1])
        print("bad temp in")
        return -np.ones(4) * np.inf
    if results_out[3][-1] < 0: 
        print("bad temp out")
        return -np.ones(4) * np.inf

    if results_in[1][-1] < 0: 
        print(results_in[1][-1])
        print("bad pressure in")
        return -np.ones(4) * np.inf
    if results_out[1][-1] < 0: 
        print("bad pressure out")
        return -np.ones(4) * np.inf

    if np.isnan(results_in).any():
        print("nan in")
        return np.ones(4) * np.inf
        
    if np.isnan(results_out).any():
        print("nan out")
        return np.ones(4) * np.inf

    # calculate differences for all parameters between inner and outer conditions, normalizes by initial guess to make all parameters on the same scale 
    difference_R = np.abs(results_in[0][-1] - results_out[0][-1])/Rs
    difference_P = np.abs(results_in[1][-1] - results_out[1][-1])/Pc 
    difference_T = np.abs(results_in[3][-1] - results_out[3][-1])/Tc
    difference_L = np.abs(results_in[2][-1] - results_out[2][-1])/Ls

    return [difference_R, difference_P, difference_L, difference_T]