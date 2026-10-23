
#Exercsie 2: Using transformation matrix to transform between local and gloobal coordinates for a composite material
#plot Q for different angles and plot strain for different stress cases and angles

#import necessary libraries and functions
import matplotlib.pyplot as plt
import numpy as np

from functions import ply, t_matrix

#number of theta values to be used for plotting
n = 31

#define arrays to store values of Q and strain for different angles
Qvalues = np.zeros((6,n)) # 6 distinct values of Q for each angle (Q11,Q22,Q12, Q16,Q26,Q66), 20 values of theta from 0 to 90 degrees
strain = np.zeros((3,n,3)) # 3 distinct values of strain for each angle (strain_x,strain_y,shear_strain), 20 values of theta from 0 to 90 degrees, 3 stress cases
strain_local = np.zeros((3,n,3)) # 3 distinct values of strain for each angle (strain_x,strain_y,shear_strain), 20 values of theta from 0 to 90 degrees, 3 stress cases
stress_local = np.zeros((3,n,3)) # 3 distinct values of stress for each angle (stress_x,stress_y,shear_stress), 20 values of theta from 0 to 90 degrees, 3 stress cases

#define theta values and stress cases

theta_list = np.linspace(-90, 90, n) #20 theta values from -90 to 90 degrees
stress_cases = [[(np.array([1,0,0])*10**6)],[(np.array([1,1,0])*10**6)],[(np.array([2, 1, 0.5])*10**6)]] #3 stress cases, 1) only one non zero, 2) two non zero stresses, 3) {2,1,0.5} stress case all in global coordinates and Mpa

#calculate local stiffness and compliance matrices for a composite material using the Halpin-Tsai equations
#S, Q_local = constitutive_matrix_HT(v_f=0.55, nu_f=0.2, E_f=76e9, E_m=4e9, nu_m=0.3, stress_state='plane_stress', psi=1) 
Q_local,S_local = ply(E1=54e9, E2=18e9, nu12=0.25, G12=9e9) #stiffness matrix for a composite material in local coordinates
#loop over theta values and compute global stiffness and compliance matrices, strains and stresses for different stress cases



for i,theta in enumerate(theta_list):
    #transform to global coordinates
    Q_global = t_matrix(theta=theta) @ Q_local @ t_matrix(theta=theta).T
    
    #save for plotting
    Qvalues[:, i] = [Q_global[0, 0], Q_global[1, 1], Q_global[0, 1], Q_global[0, 2], Q_global[1, 2], Q_global[2, 2]]
    
    #compute global strain from gloal stress using the global compliance matrix
    strain[:, i,0] = np.linalg.inv(Q_global) @  stress_cases[0][0] #apply unit stress in x-direction
    strain[:, i,1] = np.linalg.inv(Q_global) @  stress_cases[1][0] #apply unit stress in two directions (y and shear)
    strain[:, i,2] = np.linalg.inv(Q_global) @  stress_cases[2][0] #apply unit stress in shear direction

    #transform stresses and strains back to local coordinates
    strain_local[:, i, 0] = t_matrix(theta).T @ strain[:, i, 0]
    strain_local[:, i, 1] = t_matrix(theta).T @ strain[:, i, 1]
    strain_local[:, i, 2] = t_matrix(theta).T @ strain[:, i, 2]
    #local stress = T^-1 * global stress
    stress_local[:, i, 0] = np.linalg.inv(t_matrix(theta)) @ stress_cases[0][0]
    stress_local[:, i, 1] = np.linalg.inv(t_matrix(theta)) @ stress_cases[1][0]
    stress_local[:, i, 2] = np.linalg.inv(t_matrix(theta)) @ stress_cases[2][0]
    if theta == -90  or theta == 0:
        print(f"theta: {theta} degrees")
       
        print(f"strain: {strain[:, i, 2]}")
        print(f"strain_local: {strain_local[:, i, 0]}")
        print(f"stress_local: {stress_local[:, i, 0]}")
#plot 

plt.figure(figsize=(10, 6))
plt.title('Q values as function of Theta')

plt.plot(theta_list, Qvalues[0, :], label='Q11')
plt.plot(theta_list, Qvalues[1, :], label='Q22')
plt.plot(theta_list, Qvalues[2, :], label='Q12')
plt.legend()
plt.xlabel('Theta (degrees)')
plt.ylabel('Q values')

plt.figure(figsize=(10, 6))
plt.title('Q values as function of Theta')
plt.plot(theta_list, Qvalues[3, :], label='Q16')
plt.plot(theta_list, Qvalues[4, :], label='Q26')
plt.plot(theta_list, Qvalues[5, :], label='Q66')


plt.legend()
plt.xlabel('Theta (degrees)')
plt.ylabel('Q values')



#plot strains as function of theta for different stress cases
plt.figure(figsize=(10, 6))
plt.title('Global Strain as function of Theta')

plt.plot(theta_list, strain[0, :, 2], label='strain_x')
plt.plot(theta_list, strain[1, :, 2], label='strain_y')
plt.plot(theta_list, strain[2, :, 2], label='shear_strain')

plt.legend()
plt.xlabel('Theta (degrees)')
plt.ylabel('Global strain values')
#plot local strain and stress as function of theta for different stress cases
plt.figure(figsize=(10, 6))
plt.title('Local Strain values as function of Theta')


plt.plot(theta_list, strain_local[0, :, 2], label='strain_1')
plt.plot(theta_list, strain_local[1, :, 2], label='strain_2')
plt.plot(theta_list, strain_local[2, :, 2], label='shear_strain')
plt.legend()
plt.xlabel('Theta (degrees)')
plt.ylabel('Local strain values')

plt.show()




