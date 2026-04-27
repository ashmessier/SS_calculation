from load1load2 import * 
from scipy.integrate import solve_ivp
import numpy as np
from star_params import * #params0, Ms, M_intersect
from derivs import * 


def int_load1_outwards(Rs, Pc, Ls, Tc):     
    # if not (3.75 <= np.log10(Tc) <= 7.5): 
    #     return -np.inf
        
  #  M_arr = np.logspace(np.log10(M_core_init), np.log10(M_intersect), niter)[1:-1]
    params = load1(M_core_init, comp, Pc, Tc) 
    
 #   print("load1_outwards int params:", params)
    
    sol_outwards = solve_ivp(derivs, t_span = [M_core_init, M_intersect], y0 = params, method = "Radau", dense_output=True)
    return sol_outwards

def int_load2_inwards(Rs, Pc, Ls, Tc):
   # niter = int(1e5)
   # M_intersect = 0.2 * Ms

    # if not (3.75 <= np.log10(Tc) <= 7.5): 
    #     return -np.inf
        
 #   M_arr = np.logspace(np.log10(Ms), np.log10(M_intersect), niter)
   # print(M_arr[2], M_arr[-2])
    params = load2(Ms, comp, Ls, Rs) 

    sol_inwards = solve_ivp(derivs, t_span = [Ms, M_intersect], y0 = params, method = "Radau", dense_output=True)
    return sol_inwards

def difference(params):
        
    Pc, Tc, Ls, Rs = params
  #  print("step:", Pc, Tc, Ls, Rs)

    # rhoc = rho_eq(Pc, Tc, comp)
    # R = rhoc / (Tc*1e-6)**3 # calc R from T and rho 
    # logR = np.log10(R) # take log of R 
    
    # logT=np.log10(Tc)

    # if not(-8.0 < logR < 1.0):
    #     print("log R bad in difference")
    #     return np.ones(4) * np.inf
    # if not (3.75 <= logT <= 7.5): 
    #     print("log T bad")
    #     return np.ones(4) * np.inf

    results_in = int_load2_inwards(Rs, Pc, Ls, Tc).y
  #  print(results_in)
    # Ps_in = results_in[1]
    # Ts_in = results_in[3]
   # print("rho", rho_eq(Ps_in, Ts_in, comp))

    results_out = int_load1_outwards(Rs, Pc, Ls, Tc).y
    #print(results_out)

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
   # print(results_in[2][-1])
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
    
    difference_R = np.abs(results_in[0][-1] - results_out[0][-1])/Rs
    difference_P = np.abs(results_in[1][-1] - results_out[1][-1])/Pc 
    difference_T = np.abs(results_in[3][-1] - results_out[3][-1])/Tc
    difference_L = np.abs(results_in[2][-1] - results_out[2][-1])/Ls

    #print(difference_R, difference_P, difference_L, difference_T)

    return [difference_R, difference_P, difference_L, difference_T]