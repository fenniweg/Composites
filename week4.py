# Exercise 3 Compute A,B,D matrices for arbitrary layup, Ex 4 compute failures at middle of each ply using 5 different failure criteria

import numpy as np
import pandas as pd

from functions import (
    find_local_stress,
    hoffmann_criterion,
    laminate_properties,
    max_strain_criterion,
    max_stress_criterion,
    ply,
    solve_system,
    tsai_hill_criterion,
    tsai_wu_criterion,
)

###INPUTS##############################################################################
# Input the material properties for the composite material
E1 = 46e9  # Pa
E2 = 13e9  # Pa
nu12 = 0.3
G12 = 4.4e9  # Pa
t = 0.25  # thickness of each ply mm
laminate = [0, 90, 90, 0]  # laminate layup in degrees

material_strengths = [1000, 1000, 35, 150,50 ]  # [sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat_12] in Mpa


# Input the loading conditions for the laminate
M = np.array([0, 0, 5])  # moments in N
N = np.array([300, -50, 25])  # forces in N/mm

############################################################################################################

# Build the Q_matrix
Q_local = ply(E1, E2, nu12, G12)

# Build the A,B,D matrices and coordinate system for the laminates
A, B, D, z= laminate_properties(laminate, t, Q_local)
#vary the load with x 

# Solve the system for global strains and curvatures
strains, curvatures = solve_system(N, M, A, B, D)


# Initialize array to store local stresses at the bottom, middle, and top of each ply
ply_results = np.zeros((len(laminate) * 3, 4))
# Get stresses in local coordinates for bottom, middle, and top of each ply
for i, theta in enumerate(laminate):
    z_coords = [z[i], (z[i] + z[i + 1]) / 2, z[i + 1]]  # bottom, middle, and top of ply
    for j, z_coord in enumerate(z_coords):
        local_stresses = find_local_stress(strains, curvatures, z_coord, theta, Q_local)
        ply_results[i * 3 + j, 0] = z_coord
        ply_results[i * 3 + j, 1:] = local_stresses

# format the results for better readability
df = pd.DataFrame(  
    ply_results,
    columns=["z (mm)", "stress_x (MPa)", "stress_y (MPa)", "shear_stress (MPa)"],
)
    
#dict of failure criteria functions
criteria = {
    "max_stress": lambda stresses: max_stress_criterion(
        stresses, material_strengths
    ),
    "max_strain": lambda stresses: max_strain_criterion(
        stresses, material_strengths, E1, E2, G12, nu12
    ),
    "tsai_hill": lambda stresses: tsai_hill_criterion(
        stresses, material_strengths
    ),
    "tsai_wu": lambda stresses: tsai_wu_criterion(stresses, material_strengths),
    "hoffman": lambda stresses: hoffmann_criterion(stresses, material_strengths),
}
 # Initialize array to store failure results for 5 criterion (columns), for each ply 3 points(top,middle,bottom) (rows)
failure_results = np.zeros((len(laminate) * 3, 5))
#check each ply for failure using the different failure criteria
for i in range(len(laminate)):
    ply_stresses = ply_results[i * 3:(i + 1) * 3, 1:]  # get stresses for the current ply in shape of 3x3 (bottom, middle, top) three columns (stress_x, stress_y, shear_stress)
    #max_stress
    failure_results[i * 3:(i + 1) * 3, 0] = max_stress_criterion(ply_stresses,material_strengths) #result is vector containing 1 if failure occurs, 0 if not for each ply
    #max_strain
    failure_results[i * 3:(i + 1) * 3, 1] = max_strain_criterion(ply_stresses,material_strengths,E1,E2,G12,nu12) #result is vector containing 1 if failure occurs, 0 if not for each ply
    #tsai_hill
    failure_results[i * 3:(i + 1) * 3, 2] = tsai_hill_criterion(ply_stresses,material_strengths) #result is vector containing 1 if failure occurs, 0 if not for each ply
    #tsai_wu
    failure_results[i * 3:(i + 1) * 3, 3] = tsai_wu_criterion(ply_stresses,material_strengths) #result is vector containing 1 if failure occurs, 0 if not for each ply
    #hoffman
    failure_results[i * 3:(i + 1) * 3, 4] = hoffmann_criterion(ply_stresses,material_strengths) #result is vector containing 1 if failure occurs, 0 if not for each ply

df_failure = pd.DataFrame(
    failure_results,
    columns = ["max_stress", "max_strain", "tsai_hill", "tsai_wu", "hoffman"],
    index = [f"ply {i//3 + 1} {'bottom' if i%3 == 0 else 'middle' if i%3 == 1 else 'top'}" for i in range(len(laminate) * 3)]
)   

# print("Local Stresses at each ply interface (bottom, middle, top):")
# print(df)

# print("\nFailure Results for each ply interface (bottom, middle, top):")
# print(df_failure)

# print('Failure at each play middle interface (1 if failed, 0 if not):')
# print(df_failure.iloc[1::3])  # Print only the middle interface of each ply


