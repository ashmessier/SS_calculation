from opacities import *
from energies import *
from rho import *
from load1load2 import *
from derivs import *
from constants import *
from nabla_check import *
from integrate_difference import * 
from star_params import * #params0, Ms, M_intersect

import os
import numpy as np
import scipy 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp

def plot_sol(results):
    Pc, Tc, Ls, Rs = results[0]
    params = [Rs, Pc, Ls, Tc]
    int_in = int_load2_inwards(*params)
    int_out =  int_load1_outwards(*params)

    return {"m_in":int_in.t, "params_in":{"lum":int_in.y[2], "temp":int_in.y[3], "pressure":int_in.y[1], "radius":int_in.y[0]},  
            "m_out": int_out.t, "params_out":{"lum":int_out.y[2], "temp":int_out.y[3], "pressure":int_out.y[1], "radius":int_out.y[0]}}


print(f"Calculating structure for star of mass {solar_masses} solar masses")
print(f"Initial total radius guess: {params0[3]:.4g}")
print(f"Initial surface luminosity guess: {params0[2]:.4g}")
print(f"Initial core temperature guess: {params0[1]:.4g}")
print(f"Initial core pressure guess: {params0[0]:.4g}")

results = fsolve(difference, x0 = params0, full_output=True)

if results[3] != 'The solution converged.':
    raise Warning("The solution has not converged, try a different star mass?")
    print("")
   
else:
    print(f"Solution has converged after {results[1]['nfev']} iterations")    
    print("")

    plot_dict = plot_sol(results)
    
    # make plot directory 
    try:
        os.mkdir(f"plots_{solar_masses}")
        print(f"Created directory /plots_{solar_masses}")
        print("")
    except FileExistsError:
        print(f"Directory plots_{solar_masses} already exists")
        print("")

    
    fig, axes = plt.subplots(2, 2, sharex=True, dpi=750, figsize=(6,5))
    axes = axes.flatten()
    
    color_in = "teal"
    color_out = "plum"
    
    axes[0].set_title("Luminosity")
    axes[0].plot(plot_dict["m_in"], plot_dict["params_in"]["lum"], color=color_in, linewidth=2, label="Inward")
    axes[0].plot(plot_dict["m_out"], plot_dict["params_out"]["lum"], color=color_out, linewidth=2, label="Outward")
    axes[0].set_ylabel("Luminosity (ergs cm-1 s-1)")
    axes[0].axvline(0.2*Ms, color="k", linestyle='--', linewidth=0.6, label="0.2 M*")
    
    axes[1].set_title("Temperature")
    axes[1].plot(plot_dict["m_in"], plot_dict["params_in"]["temp"], color=color_in, linewidth=2)
    axes[1].plot(plot_dict["m_out"], plot_dict["params_out"]["temp"], color=color_out, linewidth=2)
    axes[1].set_ylabel("Temperature (K)")
    axes[1].axvline(0.2*Ms, color="k", linestyle='--', linewidth=0.6)
    
    axes[2].set_title("Pressure")
    axes[2].plot(plot_dict["m_in"], plot_dict["params_in"]["pressure"], color=color_in, linewidth=2)
    axes[2].plot(plot_dict["m_out"], plot_dict["params_out"]["pressure"], color=color_out, linewidth=2)
    axes[2].set_xlabel("Mass (g)")
    axes[2].set_ylabel("Pressure (dynes cm-2)")
    axes[2].axvline(0.2*Ms, color="k", linestyle='--', linewidth=0.6)    
    
    axes[3].set_title("Radius")
    axes[3].plot(plot_dict["m_in"], plot_dict["params_in"]["radius"], color=color_in, linewidth=2)
    axes[3].plot(plot_dict["m_out"], plot_dict["params_out"]["radius"], color=color_out, linewidth=2)
    axes[3].set_xlabel("Mass (g)")
    axes[3].set_ylabel("Radius (cm)")
    axes[3].axvline(0.2*Ms, color="k", linestyle='--', linewidth=0.6)
        
    axes[0].legend(frameon=True)
    fig.tight_layout()
    plt.savefig(f"plots_{solar_masses}/mass_temp_press_lum_4panel.png")

    print(f"Four panel plot saved to 'plots_{solar_masses}/mass_temp_press_lum_4panel.png'")
    print("")
    
    
    # calculate everything 
    mass_arr = np.concatenate((plot_dict["m_in"], plot_dict["m_out"][::-1]))
    lum_arr =  np.concatenate((plot_dict["params_in"]["lum"], plot_dict["params_out"]["lum"][::-1]))
    temp_arr = np.concatenate((plot_dict["params_in"]["temp"], plot_dict["params_out"]["temp"][::-1]))
    pressure_arr = np.concatenate((plot_dict["params_in"]["pressure"], plot_dict["params_out"]["pressure"][::-1]))
    radius_arr = np.concatenate((plot_dict["params_in"]["radius"], plot_dict["params_out"]["radius"][::-1]))
    
    X, Y, Z = comp
    # calc density 
    rho_arr = rho_eq(pressure_arr, temp_arr, comp)
    
    # calc energy gen rate
    energy_gen_arr = e_CNO(X, Y, Z, rho_arr, temp_arr) + e_PP(X, Y, Z, rho_arr, temp_arr)
    
    # calc radiative vs convective 
    rad_conv_arr = []
    for (m, T, P, L) in zip(mass_arr, temp_arr, pressure_arr, lum_arr):
        rad_conv_arr.append(check_nabla(comp, m, T, P, L, return_conv_rad = True))
    
    # find turning points 
    turning_masses = []
    for i, (mass, state) in enumerate(zip(mass_arr[:-1], rad_conv_arr[:-1])):
        if rad_conv_arr[i+1] != state: 
            turning_masses.append(mass)
            
    
    # calc opacity
    opacity_arr = []
    for T, rho in zip(temp_arr, rho_arr):
        logT = np.log10(T)
        logrho = np.log10(rho)
        opacity_arr.append(calc_k(logT, logrho))

    # FINISH THIS
    # calc adiabatic temp gradient 
    
    
    # calc stellar gradient dlnT/dlnP
    #stellar_gradient_arr = 
    
    # normalized 
    
    
    # plot normalized parameters 
    fig, axes = plt.subplots(1, sharex=True, dpi=750, figsize=(5,3.5))
    
    color_lum = "teal"
    color_temp = "plum"
    color_press = "darkorange"
    color_rad = "peachpuff"
    
    axes.plot(mass_arr, lum_arr/np.max(lum_arr), color=color_lum, linewidth=2, label="Luminosity")
    
    axes.plot(mass_arr, temp_arr/np.max(temp_arr), color=color_temp, linewidth=2, label="Temperature")
    
    axes.plot(mass_arr, pressure_arr/np.max(pressure_arr), color=color_press, linewidth=2, label="Pressure")
    
    axes.plot(mass_arr, radius_arr/np.max(radius_arr), color=color_rad, linewidth=2, label="Radius")
    
    # plot rad vs convective 
    axes.axvspan(0, turning_masses[2], color="coral", alpha=0.1)
    axes.axvspan(turning_masses[1], turning_masses[2], color="powderblue", alpha=0.2)
    
    #axes.semilogy()
    plt.legend(loc='upper center', 
               bbox_to_anchor=(0.5, 1.15), 
               ncol=4, fontsize=8.5)
    fig.tight_layout()
    axes.set_xlabel("Mass (g)")
    axes.set_ylabel("Normalized parameter")
    #plt.ylim(1e-3, 1e2)
    plt.savefig(f"plots_{solar_masses}/mass_temp_press_lum_one.png")
    
    print(f"Normalized parameter plot saved to 'plots_{solar_masses}/mass_temp_press_lum_one.png'")


    # calculate columns for machine readable table 
        
    mass_arr = np.concatenate((plot_dict["m_in"], plot_dict["m_out"][::-1]))
    lum_arr =  np.concatenate((plot_dict["params_in"]["lum"], plot_dict["params_out"]["lum"][::-1]))
    temp_arr = np.concatenate((plot_dict["params_in"]["temp"], plot_dict["params_out"]["temp"][::-1]))
    pressure_arr = np.concatenate((plot_dict["params_in"]["pressure"], plot_dict["params_out"]["pressure"][::-1]))
    radius_arr = np.concatenate((plot_dict["params_in"]["radius"], plot_dict["params_out"]["radius"][::-1]))
    
    X, Y, Z = comp
    
    # calc density 
    rho_arr = rho_eq(pressure_arr, temp_arr, comp)
    
    # calc energy gen rate
    energy_gen_arr = e_CNO(X, Y, Z, rho_arr, temp_arr) + e_PP(X, Y, Z, rho_arr, temp_arr)
    
    # calc radiative vs convective 
    rad_conv_arr = []
    for (m, T, P, L) in zip(mass_arr, temp_arr, pressure_arr, lum_arr):
        rad_conv_arr.append(check_nabla(comp, m, T, P, L, return_conv_rad = True))
    
    # find turning points 
    turning_masses = []
    for i, (mass, state) in enumerate(zip(mass_arr[:-1], rad_conv_arr[:-1])):
        if rad_conv_arr[i+1] != state: 
            turning_masses.append(mass)
            
    # calc opacity
    opacity_arr = []
    for T, rho in zip(temp_arr, rho_arr):
        logT = np.log10(T)
        logrho = np.log10(rho)
        opacity_arr.append(calc_k(logT, logrho))
    
    # calc adiabatic temp gradient 
    nabla_abd_arr = np.ones_like(temp_arr) * 0.4
    
    # calc stellar gradient dlnT/dlnP
    stellar_gradient_arr = []
    for (m, T, P, L) in zip(mass_arr, temp_arr, pressure_arr, lum_arr):
        stellar_gradient_arr.append(check_nabla(comp, m, T, P, L, return_conv_rad = False))
    
    #make a table of them all
    table_starparams = {"Mass":mass_arr, 
                       "Luminosity": lum_arr, 
                       "Temperature": temp_arr, 
                       "Pressure": pressure_arr, 
                       "Radius": radius_arr, 
                       "Density": rho_arr, 
                       "e": energy_gen_arr, 
                       "Energy Transport": rad_conv_arr, 
                       "Opacity": opacity_arr, 
                       "Adiabatic Temperature Gradient": nabla_abd_arr, 
                       "Stellar Gradient": stellar_gradient_arr}
    
    txt_file = f"plots_{solar_masses}/star_tab_{solar_masses}Msun.txt"
    
    # Convert to DataFrame
    df_starparams = pd.DataFrame(table_starparams)
    
    # Save to a text file
    with open(txt_file, 'w') as f:
        # to_string() keeps it as a readable table
        f.write(df_starparams.to_string(index=False))

    print(f"Saved table of star parameters to 'plots_{solar_masses}/star_tab_{solar_masses}Msun.txt'")