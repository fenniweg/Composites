'''Progressive failure analysis of Composite Laminates
Calculating FPF and LPF for different failure criteria
Printing the local stresses at the first ply failure load factor for each criterion, FPF load factor and LPF load factor for each criterion, f_critical for each ply and criterion'''

import numpy as np
import functions as func
import pandas as pd



###INPUTS##############################################################################
# Input the material properties for the composite material
E1 = 40e3  # MPa
E2 = 9.8e3  # MPa
nu12 = 0.3
G12 = 2.8e3 # MPa
t = [0.25,0.5,0.25] # 5thickness of each ply mm
laminate = [0, 45,90]  # laminate layup in degrees
load_steps = 1000

material_strengths = [1100, 600, 20, 140,70 ]  # [sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat_12] in Mpa
 # initialize failure array for each ply 0 = not failed, 1 = failed
#Loads
N = np.array([0,0,0])  #forces in N/mm [Nx, Ny, Nxy]
M = np.array([15,0, 0])  # moment around y axis [Mx, My, Mz] in N 


# Initialize array to store local stresses at the bottom, middle, and top of each ply

failure_results = np.zeros(len(laminate) * 3)  # initialize failure array for each ply 0 = not failed, 1 = failed
# Initialize results array with six rows per ply: bottom and top sigma1, sigma2, tau12.
results = np.zeros((len(laminate) * 6, 5))
fpf = np.zeros(5)  # first ply failure load factor for each criterion
lpf = np.zeros(5)  # last ply failure load factor for each criterion
f_crit = np.zeros((len(laminate)*2, 5))  # critical failure index for each ply and criterion
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


def criterion_index(stress, criterion_name):
    """Return the scalar failure index for one local stress state."""
    sigma1, sigma2, tau12 = stress
    sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12 = material_strengths
    sigma1_strength = sigma1_t if sigma1 >= 0 else sigma1_c
    sigma2_strength = sigma2_t if sigma2 >= 0 else sigma2_c

    if criterion_name == "max_stress":
        return max(
            abs(sigma1) / sigma1_strength,
            abs(sigma2) / sigma2_strength,
            abs(tau12) / tau_hat12,
        )

    if criterion_name == "max_strain":
        Q_local = func.ply(E1, E2, nu12, G12)
        epsilon1, epsilon2, gamma12 = np.linalg.solve(Q_local, stress)
        return max(
            abs(epsilon1) / (sigma1_t / E1 if epsilon1 >= 0 else sigma1_c / E1),
            abs(epsilon2) / (sigma2_t / E2 if epsilon2 >= 0 else sigma2_c / E2),
            abs(gamma12) / (tau_hat12 / G12),
        )

    if criterion_name == "tsai_hill":
        return (
            (sigma1 / sigma1_strength) ** 2
            + (sigma2 / sigma2_strength) ** 2
            - (sigma1 * sigma2) / sigma1_strength**2
            + (tau12 / tau_hat12) ** 2
        )

    if criterion_name == "tsai_wu":
        F11 = 1 / (sigma1_t * sigma1_c)
        F22 = 1 / (sigma2_t * sigma2_c)
        F12 = -0.5 * np.sqrt(F11 * F22)
        return (
            sigma1**2 / (sigma1_t * sigma1_c)
            + sigma2**2 / (sigma2_t * sigma2_c)
            + tau12**2 / tau_hat12**2
            + 2 * F12 * sigma1 * sigma2
            + sigma1 / sigma1_t
            - sigma1 / sigma1_c
            + sigma2 / sigma2_t
            - sigma2 / sigma2_c
        )

    if criterion_name == "hoffman":
        return (
            sigma1**2 / (sigma1_t * sigma1_c)
            - sigma1 * sigma2 / (sigma1_c * sigma1_t)
            + sigma2**2 / (sigma2_t * sigma2_c)
            - (sigma1_t - sigma1_c) * sigma1 / (sigma1_c * sigma1_t)
            - (sigma2_t - sigma2_c) * sigma2 / (sigma2_c * sigma2_t)
            + tau12**2 / tau_hat12**2
        )

    raise ValueError(f"Unknown failure criterion: {criterion_name}")


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
                            ply_discount[i] = 1  
                            #print(f"ply number {i + 1} failed at load factor {load_increment:.3f}") 
                            # removed from the laminate from the next step on
                            if np.sum(ply_discount) ==1: #first failure load factor
                                        #save the load factor at which the first ply fails
                                        fpf[c] = n*dM[0]
                                        
                                        #save all plie s stresses at the first failure load factor
                                        for j in range(len(laminate)):
                                            ply_stresses_j = ply_stresses[j * 3:(j + 1) * 3, :]
                                            results[j * 6:(j + 1) * 6, c] = np.concatenate(
                                                (ply_stresses_j[0, :], ply_stresses_j[2, :])
                                            )
                                            for a in[0,2]: #bottom and top surfaces
                                                f_crit[j*2+a//2, c] = criterion_index(
                                                    ply_stresses_j[a, :], criterion_name
                                                )
                                            
                                              
                                                  
        
        if np.all(ply_discount == 1):
                lpf[c] = n*dM[0]
                
                break
        if n == load_steps - 1:
            print(f"Load factor {load_increment:.2f} reached without all plies failing.")
            
df_results = pd.DataFrame(results, columns=criteria_names, index=pd.MultiIndex.from_product([["bottom", "top"], [f"ply {i + 1}" for i in range(len(laminate))], ["σ1 [MPa]", "σ2 [MPa]", "τ12 [MPa]"]]))

#add one row to save fpf and lpf for each criterion
df_results.loc[("first ply failure load factor", "", ""), :] = fpf
df_results.loc[("last ply failure load factor", "", ""), :] = lpf

#add one row to save f_critical for each ply (top and bottom) and criterion
for i in range(len(laminate)):
    df_results.loc[("critical failure index", f"ply {i + 1}", "bottom"), :] = f_crit[i * 2, :]
    df_results.loc[("critical failure index", f"ply {i + 1}", "top"), :] = f_crit[i * 2 + 1, :]

print(df_results)

#save as csv file
df_results.to_csv("progressive_failure_results.csv", index=True)

 
    
   