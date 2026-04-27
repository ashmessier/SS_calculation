from constants import a, Na, k, m_p

def rho_eq(P, T, comp):
    X, Y, Z = comp 
    mu = 4/(3+5*X-Z)

  #  mu = 4/(3 + 5*X)
    rad_pressure = (a * T**4)/3
    rho = (mu * (P-rad_pressure))/(Na * k * T)

  #  rho = (P*mu*m_p) / (k*T)
   # rho = (mu * (P-rad_pressure))/(Na * k * T)

    
    return rho 