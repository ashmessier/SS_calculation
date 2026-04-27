import numpy as np
from constants import *

def e_CNO(X, Y, Z, rho, T): 
    T9 = T/1e9
    Xcno = Z # can i do this 
    g141 = (1-2.99*T9 + 3.41*T9**2 - 2.43*T9**3)
    ecno = 8.24e25 * g141 * Xcno * X * rho * T9**(-2/3) * np.exp(-15.231*T9**(-1/3) - (T9/0.8)**2)

    return ecno

def e_PP(X, Y, Z, rho, T):
    psi = 1
    T7=T/1e7
    T9 = T/1e9
    
    if not isinstance(T, float):
        epp_arr = []
        for T7_val, T9_val, rho_val in zip(T7, T9, rho): 
            if 1.8 <= T7_val <= 2.2:  # values in book (roughly, might update later) 
                psi = 2
            elif T7_val >= 2.2:
                psi = 1.5
    
            Z1, Z2 = 1, 1
            
            f11 = np.exp(5.92e-3 *Z1*Z2*(1*rho_val/T7_val**3)**1/2)
            
            g11 = (1 + 3.82*T9_val + 1.51*T9_val**2 + 0.144*T9_val**3 - 0.0114*T9_val**4)
            epp = 2.57e4 * psi*f11*g11*rho_val*X**2*T9_val**(-2/3) * np.exp(-3.381/(T9_val**(1/3)))
            epp_arr.append(epp)
            return epp_arr

    else: 
        if 1.8 <= T7 <= 2.2:  # values in book (roughly, might update later) 
            psi = 2
        elif T7 >= 2.2:
            psi = 1.5
    
        Z1, Z2 = 1, 1
        
        f11 = np.exp(5.92e-3 *Z1*Z2*(1*rho/T7**3)**1/2)
        
        g11 = (1 + 3.82*T9 + 1.51*T9**2 + 0.144*T9**3 - 0.0114*T9**4)
        epp = 2.57e4 * psi*f11*g11*rho*X**2*T9**(-2/3) * np.exp(-3.381/(T9**(1/3)))

        return epp 
