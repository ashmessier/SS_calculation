from constants import Rsun, Msun, Lsun 

solar_masses = input("Input desired star mass in solar masses: ")

Ms = Msun * float(solar_masses)
M_core_init = 1e10
M_intersect = 0.2 * Ms
niter = int(1e3)

Rs_guess = Rsun*(float(solar_masses)**(0.75)) #* 1.67 # do scaling relations 
Pc_guess = 1.5e17 
Ls_guess = Lsun * (float(solar_masses)**(3.5))  #* 15.33 # scaling relations 
Tc_guess = 2e7 

comp = (0.70, 0.28, 0.02)

params0 = [Pc_guess, Tc_guess, Ls_guess, Rs_guess]