from constants import a, Na, k, m_p

def rho_eq(P, T, comp):

    """
    Calulate density given pressure, temperature, and composition using the equation of state from radiation pressure and the ideal gas law.
    All parameters in their respective cgs units. 

    inputs: 
    P (float): Pressure at certain mass  
    T (float): Temperature at certain mass 
    comp: Tuple(Float X, Float Y, Float Z): fractional composition of the star in (fracction hydrogen, fraction helium, fraction metals) 

    returns: 
    rho (float): Density of the star 
    """
    
    X, Y, Z = comp 
    mu = 4/(3+5*X-Z)

    rad_pressure = (a * T**4)/3
    rho = (mu * (P-rad_pressure))/(Na * k * T)

    
    return rho 