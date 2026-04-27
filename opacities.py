import numpy as np
import scipy 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from constants import *


# opacities
star_comp_table = "star_comp_key.txt" # get edited table w indices 

star_comp_csv = pd.read_csv(star_comp_table, delimiter=" ") # read in 

def get_table_number(X, Y, Z): 
    X_filt = star_comp_csv[star_comp_csv['x_val'] == "X="+X] 
    Y_filt = X_filt[X_filt['y_val'] == "Y="+Y]
    
    if np.shape(Y_filt)[0] != 1:  # if multiple Z vals left; use specified Z value, 
        Z_filt = Y_filt[Y_filt['z_val'] == "Z="+Z]
        return int(Z_filt["number"].iloc[0])
    else: 
        return int(Y_filt["number"].iloc[0]) # else just return the vals filtered on X and Y 
        
table_indx = get_table_number('0.7000', '0.2800', '0.0000') # specified vals from class 

tables_file = '/Users/asmessier/Desktop/JHU/Spring2026/stars/stellar_structure_calc/opacites_hw2.txt'

tab_indxs = {} # initialize dict for table row lines 

with open(tables_file, 'r') as file: # go through all lines in big csv 
    table, tab_start_val, tab_end_val = [0, 0, 0]
    for i, line in enumerate(file): 
        if "TABLE" in line:
            table = line.split("#")[1].split("     ")[0] # gives number of table 
        if "logT" in line:  # gives first row of table 
            tab_start_val = i
        if "8.70" in line: 
            tab_end_val = i # gives last row of table 

        tab_indxs[int(table)] = (tab_start_val, tab_end_val)
        
def get_table(table_indx):
    final_k_table = pd.read_csv(tables_file, 
        skiprows=tab_indxs[table_indx][0], 
        skipfooter=int(9576-tab_indxs[table_indx][1]-1),
        sep='\s+', engine='python') # opens final table skipping necessary rows calculated prior 
    return final_k_table

final_k_table = get_table(table_indx)

# using RegularGridInterpolator function to interpolate table -> need list of x and y coords of table 
logT_coords = np.array(final_k_table["logT"])  # get logT (first column) 
logR_coords = [float(val) for val in final_k_table.columns.tolist()[1:]] # gets list of first row of table 

final_vals = final_k_table.drop(columns=final_k_table.columns[:1]) # do this to get header 

new_header = final_vals.iloc[1] #grab the first column of final vals 
final_vals = final_vals[0:] #take the entire data minus the header row
final_vals.columns = new_header #set the header row as the df header
vals_arr = final_vals.to_numpy() # convert to array 

interpolator = RegularGridInterpolator((logT_coords, 
                        logR_coords), vals_arr, method='linear') # make interpolator object 

def calc_k(logT, logrho): # use relation in doccument 
    """
    The logarithm of the Rosseland mean opacity [cm**2/g] as a function
    of log(T) for columns of constant log(R), where
 
    R=density[g/cm**3]/T6**3, T6=1.e-6*T[degrees]
    log(T) range: 70 values from 3.75 to 8.70
    log(R) range: 19 values from -8.0 to +1.0
    
    """
    
    T = 10**logT # linea/r T 
    rho = 10**logrho # linear rho
    R = rho / (T*1e-6)**3 # calc R from T and rho 
    logR = np.log10(R) # take log of R 
    
  #  print(f"Calculating opacity for logR = {logR:.2f}; logT = {logT:.2f}")

    # make sure vals for T and R are in range of table 
    if not (3.75 <= logT <= 7.5): 
        raise Warning("Log(T) out of bounds, adjust T")

    if not (-8.0 <= logR <= 1): 
        raise Warning("Log(R) out of bounds, adjust rho")
     
    interpolated_k = interpolator([logT, logR]) 
   # print(float(10**interpolated_k[0]))
    return float(10**interpolated_k[0])
    
