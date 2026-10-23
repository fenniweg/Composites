'''Progressive failure analysis of Composite Laminates'''

import numpy as np
import functions as func
import pandas as pd
import sys


###INPUTS##############################################################################
# Input the material properties for the composite material
E1 = 40e3  # MPa
E2 = 9.8e3  # MPa
nu12 = 0.3
G12 = 2.8e3 # MPa
t = [0.25,0.5,0.25] # 5thickness of each ply mm
laminate = [0, 45,90]  # laminate layup in degrees
load_steps = 10000

material_strengths = [1100, 600, 20, 140,70 ]  # [sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat_12] in Mpa
 # initialize failure array for each ply 0 = not failed, 1 = failed
#Loads
N = np.array([0,0,0])  #forces in N/mm [Nx, Ny, Nxy]
M = np.array([15,0, 0])  # moment around y axis [Mx, My, Mz] in N 


# Initialize array to store local stresses at the bottom, middle, and top of each ply

failure_results = np.zeros(len(laminate) * 3)  # initialize failure array for each ply 0 = not failed, 1 = failed
#Initialize results array to store reulst for all criteria failure (columns) for all plies stresses of failure (rows) sigma1,sigma2,tau12   
results = np.zeros((len(laminate) *4*2+2, 5))  # initialize results array to store results for all criteria failure (columns) 
#rows for each lamina top and bottom, sigma1,ssigma2,tau12, critical load, plus save first failure load and last failure load load value
   #dict of failure criteria functions
criteria = {
    "max_stress": lambda stresses: func.max_stress_criterion(
        stresses, material_strengths
    ),
    "max_strain": lambda stresses: func.max_strain_criterion(
        stresses, material_strengths, E1, E2, G12, nu12
    ),
    "tsai_hill": lambda stresses: func.tsai_hill_criterion(
        stresses, material_strengths
    ),
    "tsai_wu": lambda stresses: func.tsai_wu_criterion(stresses, material_strengths),
    "hoffman": lambda stresses: func.hoffmann_criterion(stresses, material_strengths),
}
criteria_names = list(criteria.keys())
dN = N/load_steps
dM = M/load_steps
#looping over stress criteria
for c,criterion_name in enumerate(criteria_names):
    print(f"\nEvaluating failure criterion: {criterion_name}")
    ply_results = np.zeros((len(laminate) * 3, 4)) # initialize failure array for each ply 0 = not failed, 1 = faile
    ply_stresses = np.zeros((len(laminate) * 3, 3)) # Initialize array to store local stresses at the bottom, middle, and top of each ply

    Q_local = func.ply(E1, E2, nu12, G12) #calculate the local stiffness matrix for a single ply

    ply_discount = np.zeros(len(laminate)) #ply discount array to keep track of failed plies, 0 = not failed, 1 = failed
    strains_total = np.zeros(3) #initialize total strains to zero
    curvatures_total = np.zeros(3) #initialize total curvatures to zero
    #loop over load steps to incrementally increase the load and check for ply failure
    for n in range(load_steps):  # loop over load steps to incrementally increase the load and check for ply failure
        
        #current load factor
        load_increment = (n + 1) / load_steps   
    
        #get the A,B,D matrices and z coordinates for each ply interface with failed plies discounted
        A,B,D,z = func.laminate_properties_discount(laminate,t, Q_local, ply_discount) #matrices for laminate with failed plies fully  discounted
        
        #solve the system for global strains and curvatures
        dstrains, dcurvatures = func.solve_system(dN, dM, A, B, D)
        strains_total += dstrains
        curvatures_total += dcurvatures

        strains,curvatures = strains_total, curvatures_total  
          # Exit the program after printing the strains for debugging

        #solve for local stresses at the bottom, middle, and top of each ply
        for i, theta in enumerate(laminate):
            z_coords = [z[i], (z[i] + z[i + 1]) / 2, z[i + 1]]  # bottom, middle, and top of ply
            for j, z_coord in enumerate(z_coords):
                local_stresses = func.find_local_stress(strains, curvatures, z_coord, theta, Q_local)
                ply_results[i * 3 + j, 0] = z_coord
                ply_stresses[i * 3 + j, :] = local_stresses
            
            
        for i in range(len(laminate)):
            if ply_discount[i] == 1:
                continue  # Skip this ply if it has failed
            ply_stresses_i = ply_stresses[i * 3:(i + 1) * 3, :]  # get stresses for the current ply in shape of 3x3 (bottom, middle, top) three columns (stress_x, stress_y, shear_stress)
            surface_failure = criteria[criterion_name](ply_stresses_i)[[0, 2]]

            # A failed ply is fully degraded; only its bottom and top surfaces
            # are used to trigger failure.
            if np.any(surface_failure):
                            ply_discount[i] = 1   # removed from the laminate from the next step on
                            if np.sum(ply_discount) ==1: #first failure load factor
                                        m_min = n*dM[0]
                                        print(f"First ply failure at load M = {m_min:.3f} N-mm/mm")
                            #print(f"Ply {i + 1} failed at load factor {load_increment:.3f}")
                            results[i * 6:(i + 1) * 6, c] = np.concatenate(
                                (ply_stresses_i[0, :], ply_stresses_i[2, :])
                            )
        
        if np.all(ply_discount == 1):
                m_max = n*dM[0]
                print(f"Last ply failure at load M = {m_max:.3f} N-mm/mm")
                break
        if n == load_steps - 1:
            print(f"Load factor {load_increment:.2f} reached without all plies failing.")
            
df_results = pd.DataFrame(results, columns=criteria_names)
print(df_results)

   
    
   