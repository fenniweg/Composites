
#Functions.py
#Functions used for calculations of fiber composites 

import numpy as np
import sympy as sp


def solve_system(n,m, A,B,D):
    """
    Solves the system of equations for a laminate given the applied forces n and moments m  and the A,B,D matrices and returns two arrays for strains and curvatures.
    """
    
    #left side n,m 1x6
    left_side = np.concatenate((n, m), axis=0)

    #matrix A,B;B,D 6x6
    matrix = np.block([[A, B], [B, D]])

    
    # Solve the system of equations
    solution = np.linalg.solve(matrix, left_side)
    strains = solution[:3] #strains in x,y,shear directions
    curvatures = solution[3:6] #curvatures in x,y,shear directions
    
    return strains,curvatures



def t_matrix(theta):
    """
    Returns the transformation matrix for a given angle theta (in degrees).
    """
    theta_rad = np.radians(theta)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)

    T = np.array([[c**2, s**2, -2*c*s],
                  [s**2, c**2, 2*c*s],
                  [s*c, -c*s, c**2 - s**2]])
    
    return T
def constitutive_matrix(v_f,nu_f, E_f, E_m,nu_m,stress_state):
    """
    Returns the constitutive matrix S and Q for a composite material given its fiber volume fraction, fiber and matrix properties,
    and stress state.
    """

    #Young's modulus of the composite material
    E1 =  (1-v_f)*E_m+v_f*E_f
    E2 = (v_f/E_f+(1-v_f)/E_m)**-1
    E3 = E2

    
        

    #Bulk modulus of the composite material
    Kf = E_f/(3*(1-2*nu_f))
    Km = E_m/(3*(1-2*nu_m))
    K = (v_f/Kf + ((1-v_f)/Km))**-1

    

    #Poisson's ratio of the composite material
    nu12 = nu_f*v_f + nu_m*(1-v_f)
    nu21 = (nu_f*v_f+(1-v_f)*nu_m)*E2/E1
    nu23 = 1-nu21-(E2/(3*K))

    #Shear modulus of the composite material
    G23 = E2/(2*(1+nu23))
    Gm = E_m/(2*(1+nu_m))
    Gf = E_f/(2*(1+nu_f))
    G12 = (v_f/Gf+(1-v_f)/Gm)**-1

    #Poisson's ratio of the composite material
    nu21= E2*nu12/E1
    
    nu13 = nu12

    G13 = G12

    # print("E1:", E1)
    # print("E2:", E2)
    # print("E3:", E3)
    # print("nu12:", nu12)
    # print("nu13:", nu13)
    # print("nu21:", nu21)
    # print("nu23:", nu23)
    # print("nu31:", nu31)
    # print("nu32:", nu23)
    # print("G12:", G12)
    # print("G23:", G23)
    # print("G13:", G13)
      


    #Assemble the constitutive matrix S and Q for plane stress and plane strain conditions,  matrix S for 3D

    if stress_state == "plane_stress":
        matrix_S = np.array([[1/E1, -nu12/E2, 0],  
                         [-nu12/E1, 1/E2, 0],
                         [0, 0, 1/G12]])

        matrix_Q = np.linalg.inv(matrix_S)

    else:

        matrix_S = np.array([[1/E1, -nu12/E2, -nu13/E3, 0, 0, 0],  
                         [-nu12/E1, 1/E2, -nu23/E3, 0, 0, 0],
                         [-nu13/E1, -nu23/E2, 1/E3, 0, 0, 0],
                            [0, 0, 0, 1/G23, 0, 0],
                            [0, 0, 0, 0, 1/G13, 0],
                            [0, 0, 0, 0, 0, 1/G12]])
        matrix_Q = 0



    return matrix_S,matrix_Q


def constitutive_matrix_HT(v_f,nu_f, E_f, E_m,nu_m,stress_state,psi):
    """
    Returns the constitutive matrix S and Q for a composite material given its fiber volume fraction, fiber and matrix properties,
    stress state, and Halpin-Tsai parameter psi.
    """
    #Young's modulus of the composite material
    E1 =  (1-v_f)*E_m+v_f*E_f
    E_f2 = E_f

    eta1=(E_f2/E_m-1)/(E_f2/E_m+psi)

    E2 = E_m*(1+psi*v_f*eta1)/(1-eta1*v_f)
    E3 = E2

    #Bulk modulus of the composite material
    Kf = E_f/(3*(1-2*nu_f))
    Km = E_m/(3*(1-2*nu_m))
    K = (v_f/Kf + ((1-v_f)/Km))**-1
    
        
    
    #Poisson's ratio of the composite material
    nu12 = nu_f*v_f + nu_m*(1-v_f)
    nu21 = (nu_f*v_f+(1-v_f)*nu_m)*E2/E1
    nu23 = 1-nu21-(E2/(3*K))
    

    
        

    #Bulk modulus of the composite material
    Kf = E_f/(3*(1-2*nu_f))
    Km = E_m/(3*(1-2*nu_m))
    K = (v_f/Kf + ((1-v_f)/Km))**-1

    

    #Poisson's ratio of the composite material
    nu12 = nu_f*v_f + nu_m*(1-v_f)
    nu21 = (nu_f*v_f+(1-v_f)*nu_m)*E2/E1
    nu23 = 1-nu21-(E2/(3*K))

     #Shear modulus of the composite material
    G23 = E2/(2*(1+nu23))
    Gm = E_m/(2*(1+nu_m))
    Gf = E_f/(2*(1+nu_f))
    psi2 = psi

    eta2 = (Gf/Gm-1)/(Gf/Gm+psi2)
    G12 = Gm*(1+psi*eta2*v_f)/(1-eta2*v_f)
    
        #Poisson's ratio of the composite material
    nu21= E2*nu12/E1
        
    nu13 = nu12
    
    G13 = G12
    
    # print("E1:", E1)
    # print("E2:", E2)
    # print("E3:", E3)
    # print("nu12:", nu12)
    # print("nu13:", nu13)
    # print("nu21:", nu21)
    # print("nu23:", nu23)
    # print("nu31:", nu31)
    # print("nu32:", nu23)
    # print("G12:", G12)
    # print("G23:", G23)
    # print("G13:", G13)

    #AssE_mble the constitutive matrix S and Q for plane stress and plane strain conditions,  matrix S for 3D

    if stress_state == "plane_stress":
        matrix_S = np.array([[1/E1, -nu21/E2, 0],  
                         [-nu12/E1, 1/E2, 0],
                         [0, 0, 1/G12]])

        matrix_Q = np.linalg.inv(matrix_S)

    else:

        matrix_S = np.array([[1/E1, -nu12/E2, -nu13/E3, 0, 0, 0],  
                            [-nu12/E1, 1/E2, -nu23/E3, 0, 0, 0],
                            [-nu13/E1, -nu23/E2, 1/E3, 0, 0, 0],
                            [0, 0, 0, 1/G23, 0, 0],
                            [0, 0, 0, 0, 1/G13, 0],
                            [0, 0, 0, 0, 0, 1/G12]])
        matrix_Q = 0


    return matrix_S,matrix_Q


def ply(E1,E2,nu12,G12):

    """
    Returns the inverse stiffness matrix Q for a ply given its material properties.
    """
    nu21 = (E2 * nu12) / E1
    S11 = 1/E1
    S12 = -nu21/E2
    S21 = -nu12/E1
    S22 = 1/E2
    S66 = 1/G12

    S_local = np.array([[S11, S12, 0],
                        [S21, S22, 0],
                        [0, 0, S66]])
    Q_local = np.linalg.inv(S_local)
    
    return Q_local


def laminate_properties(laminate, t, Q_local):
    """
    
    Gives the coordinate system with z=0 at the mid-plane of the laminate. Returns A, B, D matrices and z coordinates for each ply interface.
    """
    #check if t is a scalar or an array, if scalar make it an array of the same length as laminate
    if np.isscalar(t):
        t = np.full(len(laminate), t)
    
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    z = np.zeros(len(laminate)+1) #initialize array to store z coordinates for each ply interface
    
    full_thickness = sum(t)
    
        
    z[0] = -full_thickness/2.0 #bottom of laminate
    for k in range(len(laminate)):
            z[k+1] = z[k] + t[k] #top of each ply
    
    for i,theta in enumerate(laminate):
        
        #compute z coordinates for top and bottom of each ply
            z_bottom = -full_thickness/2 + sum(t[:i])
            z_top = z_bottom + t[i]
            
            #transform to global coordinates
            Q_global = t_matrix(theta=theta) @ Q_local @ t_matrix(theta=theta).T
        
            #compute A,B,D matrices
            A += Q_global * (z_top - z_bottom)
            B += 0.5 * Q_global * (z_top**2 - z_bottom**2)
            D += (1/3) * Q_global * (z_top**3 - z_bottom**3)
        #check each element of A,B,D matrices and set to zero if less than a threshold value
    for element in [A, B, D]:
            element[np.abs(element) < 1e-5] = 0
    
         
        
    return A, B, D, z


def find_local_stress(strains,curvatures,z_coord,theta,Q_local):
    """
    Computes the local strains for each ply given the global strains and curvatures.
    """    
    global_strains = strains + z_coord * curvatures
    local_strains = t_matrix(theta=theta).T @ global_strains
    local_stresses = Q_local @ local_strains
   
    return local_stresses



def max_stress_criterion(local_stresses,material_strengths):
    """Uses the maximum stress criterion to determine if the ply has failed based on the local stresses (sigma1, sigma2, tau12) at top,middle,bottom based on
    material strenghts (sigmah1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12) for the ply and return a vector of 1s and 0s for each ply interface, 1 if failed, 0 if not failed.
    """
    failure = np.zeros(len(local_stresses))
    sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12 = material_strengths
    for i in range(len(local_stresses)):
        sigma1, sigma2, tau12 = local_stresses[i]
        if (sigma1 > sigma1_t) or (sigma1 < -sigma1_c) or (sigma2 > sigma2_t) or (sigma2 < -sigma2_c) or (np.abs(tau12 )> tau_hat12):
            failure[i] =1

    return failure

def max_strain_criterion(local_stresses,material_strengths,E1,E2,G12,nu12):
    """Uses the maximum strain criterion to determine if the ply has failed based on the local strains (epsilon1, epsilon2, gamma12) at top,middle,bottom based on
    material strains (epsilon1_t, epsilon1_c, epsilon2_t, epsilon2_c, gamma_hat12) for the ply and return a vector of 1s and 0s for each ply interface, 1 if failed, 0 if not failed.
    """
    failure = np.zeros(len(local_stresses))
    material_strains = [material_strengths[0]/E1, material_strengths[1]/E1, material_strengths[2]/E2, material_strengths[3]/E2, material_strengths[4]/G12]
    Q_local = ply(E1, E2, nu12, G12)
    local_strains = local_stresses @ np.linalg.inv(Q_local)  # Convert local stresses to local strains using the inverse of the stiffness matrix
    epsilon1_t, epsilon1_c, epsilon2_t, epsilon2_c, gamma_hat12 = material_strains
    
    for i in range(len(local_strains)):
        epsilon1, epsilon2, gamma12 = local_strains[i]
        if (epsilon1 > epsilon1_t) or (epsilon1 < -epsilon1_c) or (epsilon2 > epsilon2_t) or (epsilon2 < -epsilon2_c) or (np.abs(gamma12 )> gamma_hat12):
            failure[i] =1

    return failure

def tsai_hill_criterion(local_stresses, material_strengths):
    """Uses the Tsai-Hill criterion to determine if the ply has failed based on the local stresses (sigma1, sigma2, tau12) at top,middle,bottom based on
    material strengths (sigmah1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12) for the ply and return a vector of 1s and 0s for each ply interface, 1 if failed, 0 if not failed.
    """
    failure = np.zeros(len(local_stresses))
    sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12 = material_strengths
    for i in range(len(local_stresses)):
        sigma1, sigma2, tau12 = local_stresses[i]
        #check if compression or tension failure based on the sign of the stress to choose correct strength values
        if sigma1 > 0:
            sigma1_material = sigma1_t
        else:
            sigma1_material = sigma1_c
        if sigma2 > 0:
            sigma2_material = sigma2_t
        else:
            sigma2_material = sigma2_c

        
        #Tsai-Hill criterion
        if (sigma1/sigma1_material)**2 + (sigma2/sigma2_material)**2 - (sigma1*sigma2)/(sigma1_material**2) + (tau12/tau_hat12)**2 > 1:
            failure[i] = 1

    return failure


def tsai_wu_criterion(local_stresses, material_strengths):
    """Uses the Tsai-Wu criterion to determine if the ply has failed based on the local stresses (sigma1, sigma2, tau12) at top,middle,bottom based on
    material strengths (sigmah1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12) for the ply and return a vector of 1s and 0s for each ply interface, 1 if failed, 0 if not failed.
    """
    failure = np.zeros(len(local_stresses))
    sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12 = material_strengths
    F12_star = -0.5 
    for i in range(len(local_stresses)):
        sigma1, sigma2, tau12 = local_stresses[i]
        F11 = 1/(sigma1_t*sigma1_c)
        F22 = 1/(sigma2_t*sigma2_c)
        F12 = F12_star*np.sqrt(F11*F22)
        if sigma1**2/(sigma1_t*sigma1_c) + sigma2**2/(sigma2_t*sigma2_c) + tau12**2/(tau_hat12**2) + 2*F12*(sigma1*sigma2) + (sigma1/sigma1_t) - (sigma1/sigma1_c) + (sigma2/sigma2_t) - (sigma2/sigma2_c) > 1:
            failure[i] = 1
            failure_index = i

    return failure

def hoffmann_criterion(local_stresses, material_strengths):
    """Uses the Hoffmann criterion to determine if the ply has failed based on the local stresses (sigma1, sigma2, tau12) at top,middle,bottom based on
    material strengths (sigmah1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12) for the ply and return a vector of 1s and 0s for each ply interface, 1 if failed, 0 if not failed.
    """
    failure = np.zeros(len(local_stresses))
    sigma1_t, sigma1_c, sigma2_t, sigma2_c, tau_hat12 = material_strengths
    for i in range(len(local_stresses)):
        sigma1, sigma2, tau12 = local_stresses[i]
        if sigma1**2/(sigma1_t*sigma1_c) -sigma1*sigma2/(sigma1_c*sigma1_t) + sigma2**2/(sigma2_t*sigma2_c) - (sigma1_t-sigma1_c)*sigma1/(sigma1_c*sigma1_t)-(sigma2_t-sigma2_c)*sigma2/(sigma2_c*sigma2_t)+tau12**2/tau_hat12**2 > 1:
            failure[i] = 1

    return failure





def laminate_properties_discount(laminate, t, Q_local,ply_discount):
    """
    t is an array of thickness for each ply
    Gives the coordinate system with z=0 at the mid-plane of the laminate. Includes ply discount for failed laminates Returns A, B, D matrices and z coordinates for each ply interface.
    failure is a an array containing 1 for failed plies and 0 for non-failed plies (len(laminate)). Failed plies are discounted from the A,B,D matrices.
    """
    #check if t is a scalar or an array, if scalar make it an array of the same length as laminate
    if np.isscalar(t):
        t = np.full(len(laminate), t)
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    z = np.zeros(len(laminate)+1) #initialize array to store z coordinates for each ply interface

    full_thickness = sum(t)

    
    z[0] = -full_thickness/2.0 #bottom of laminate
    for k in range(len(laminate)):
        z[k+1] = z[k] + t[k] #top of each ply

    for i,theta in enumerate(laminate):
        if (ply_discount[i] == 1):  # Check if any the ply has failed and is discounted
            continue  # Skip this ply if it has failed
        #compute z coordinates for top and bottom of each ply
        z_bottom = -full_thickness/2 + sum(t[:i])
        z_top = z_bottom + t[i]
        
        #transform to global coordinates
        Q_global = t_matrix(theta=theta) @ Q_local @ t_matrix(theta=theta).T
    
        #compute A,B,D matrices
        A += Q_global * (z_top - z_bottom)
        B += 0.5 * Q_global * (z_top**2 - z_bottom**2)
        D += (1/3) * Q_global * (z_top**3 - z_bottom**3)
    #check each element of A,B,D matrices and set to zero if less than a threshold value
    for element in [A, B, D]:
        element[np.abs(element) < 1e-5] = 0

     
    
    return A, B, D, z



def deflection_plate(a, b, q, D,n_fourier,m_fourier,alpha1,beta1,alpha2,beta2):

    """ 
    Calculates the deflection at the center of a simply supported rectangular plate under uniform load using classical plate theory.
    a: length of the plate in x direction (m),b: length of the plate in y direction (m), p0: uniform load (Pa), D: bending stiffness matrix, 
    n, m : number of terms
    NOT USED ANYMORE, USE orthotropic_plate_arb() INSTEAD
    """

    D11 = D[0, 0]
    D12 = D[0, 1]
    D66 = D[2, 2]
    D22 = D[1, 1]

    print(f"D11: {D11}, D12: {D12}, D22: {D22}, D66: {D66}")

    #load at middle of plate 
    x = a / 2
    y  = b / 2

    
    w = 0
    #fourier 
    for n in range(1, n_fourier+1):

        for m in range(1, m_fourier+1):            
            w += np.sin(m*np.pi*x/a)*np.sin(n*np.pi*y/b)/(m*n*(D11*(m*b/a)**4+2*(D12+2*D66)*(m*n*b/a)**2+D22*n**4))
            
    w0 = 16 *q*b**4/np.pi**6
    w = w0 * w
    
    return w    


def moments_plate(a, b, q, D,n_fourier,m_fourier):
    """ calculate moments Mx,My and Mxy  on a simply supported rectangular plate under uniform load using classical plate theory.
    a: length of the plate in x direction (m),b: length of the plate in y direction (m), p0: uniform load (Pa), D: bending stiffness matrix, n, m : number of terms
    Mx and My are evaluated at the center of the plate (a/2, b/2) and Mxy is evaluated at the corners of the plate (a, 0)
    NOT USED ANYMORE, USE orthotropic_plate_arb() INSTEAD"""


    D11 = D[0, 0]
    D12 = D[0, 1]
    D16 = D[0, 2]
    D26 = D[1, 2]
    D22 = D[1, 1]
    D66 = D[2, 2]

    #find double derivative of deflection w with respect to x and y using sympy
    x, y = sp.symbols('x y')

    w = 0
    for n in range(1, n_fourier+1):
    
            for m in range(1, m_fourier+1):            
                w += sp.sin(m*sp.pi*x/a)*sp.sin(n*sp.pi*y/b)/(m*n*(D11*(m*b/a)**4+2*(D12+2*D66)*(m*n*b/a)**2+D22*n**4))
                
    w0 = 16 *q*b**4/sp.pi**6
    w = w0 * w    

    wxx = sp.diff(sp.diff(w, x), x)
    wyy = sp.diff(sp.diff(w, y), y)
    wxy = sp.diff(sp.diff(w, x), y)      

    M_x = -(D11 * wxx + D12 * wyy + 2* D16 * wxy)
    M_x = M_x.subs({x: a/2, y: b/2})
    M_x = M_x.evalf()

    M_y = -(D12 * wxx+ D22 * wyy + 2* D26 * wxy)
    M_y = M_y.subs({x: a/2, y: b/2})
    M_y = M_y.evalf()

    #the absolute value of Mxy is evaluated at the corners
    M_xy_abs = 0
    M_xy_func = -(D16 * wxx+ D26 * wyy + 2* D66 * wxy)
    for x_sub in [0, a]:
        for y_sub in [0, b]:
            M_xy = M_xy_func.subs({x: x_sub, y: y_sub})
            M_xy = M_xy.evalf()
            M_xy_abs  = max(M_xy_abs, np.abs(M_xy))
            

      
    return M_x,M_y,M_xy_abs

    
def orthropic_plate(a, b, q, D,n_fourier,m_fourier):

    """ 
    Calculates the deflection at the center of a simply supported rectangular plate under uniform load using classical plate theory.
    Calculates Moments Mx, My at the center of the plate and Mxy at the corners of the plate returing the maximum absolute value under unidorm loading.
    a: length of the plate in x direction (m),b: length of the plate in y direction (m), p0: uniform load (Pa), D: bending stiffness matrix, 
    n, m : number of terms
    USE orthotropic_plate_arb() INSTEAD FOR ARBITRARY LOADS 
    """

    D11 = D[0, 0]
    D12 = D[0, 1]
    D66 = D[2, 2]
    D22 = D[1, 1]
    D16 = D[0, 2]
    D26 = D[1, 2]

    print(f"D11: {D11}, D12: {D12}, D22: {D22}, D66: {D66}")

    x, y = sp.symbols('x y')
    
    w = 0
    for n in range(1, n_fourier+1):
        
                for m in range(1, m_fourier+1):            
                    w += sp.sin(m*sp.pi*x/a)*sp.sin(n*sp.pi*y/b)/(m*n*(D11*(m*b/a)**4+2*(D12+2*D66)*(m*n*b/a)**2+D22*n**4))
                    
    w0 = 16 *q*b**4/sp.pi**6
    w = w0 * w 

    deflection = w.subs({x:a/2, y:b/2}).evalf()

    wxx = sp.diff(sp.diff(w, x), x)
    wyy = sp.diff(sp.diff(w, y), y)
    wxy = sp.diff(sp.diff(w, x), y)      
    
    M_x = -(D11 * wxx + D12 * wyy + 2* D16 * wxy)
    M_x = M_x.subs({x: a/2, y: b/2}).evalf()
    
    
    M_y = -(D12 * wxx+ D22 * wyy + 2* D26 * wxy)
    M_y = M_y.subs({x: a/2, y: b/2}).evalf()
    
        #the absolute value of Mxy is evaluated at the corners
    M_xy_abs = 0
    M_xy_func = -(D16 * wxx+ D26 * wyy + 2* D66 * wxy)
    for x_sub in [0, a]:
            for y_sub in [0, b]:
                M_xy = M_xy_func.subs({x: x_sub, y: y_sub}).evalf()
                M_xy_abs  = max(M_xy_abs, np.abs(M_xy))
                
    
          
    

    
    return deflection, M_x, M_y, M_xy_abs



def orthropic_plate_arb(a, b, q, D,n_fourier,m_fourier,alpha1,beta1,alpha2,beta2):

    """ 
    Calculates the deflection at the center of a simply supported rectangular plate under arbitrary load using classical plate theory. 
    Calculates Moments Mx, My at the center of the plate and Mxy at the corners of the plate returing the maximum absolute value.
    a: length of the plate in x direction (m),b: length of the plate in y direction (m), p0: uniform load (Pa), D: bending stiffness matrix, alpha1, beta1, alpha2, beta2: load distribution parameters, n, m : number of terms
    """

    D11 = D[0, 0]
    D12 = D[0, 1]
    D66 = D[2, 2]
    D22 = D[1, 1]
    D16 = D[0, 2]
    D26 = D[1, 2]

    x, y = sp.symbols('x y')
    
    w = 0
    for n in range(1, n_fourier+1):
                for m in range(1, m_fourier+1):  
                    qmn = 4 *q/(sp.pi**2*m*n)*(sp.cos(m*sp.pi*alpha2)-sp.cos(m*sp.pi*alpha1))*(sp.cos(n*sp.pi*beta2)-sp.cos(n*sp.pi*beta1))
                    w += qmn*sp.sin(m*sp.pi*x/a)*sp.sin(n*sp.pi*y/b)/(D11*(m*sp.pi/a)**4+2*(D12+2*D66)*(m*sp.pi/a)**2*(n*sp.pi/b)**2+D22*(n*sp.pi/b)**4)
                    
    

    deflection = w.subs({x:a/2, y:b/2}).evalf()

    wxx = sp.diff(sp.diff(w, x), x)
    wyy = sp.diff(sp.diff(w, y), y)
    wxy = sp.diff(sp.diff(w, x), y)      
    
    M_x = -(D11 * wxx + D12 * wyy + 2* D16 * wxy)
    M_x = M_x.subs({x: a/2, y: b/2}).evalf()
    
    
    M_y = -(D12 * wxx+ D22 * wyy + 2* D26 * wxy)
    M_y = M_y.subs({x: a/2, y: b/2}).evalf()
    
        #the absolute value of Mxy is evaluated at the corners
    M_xy_abs = 0
    M_xy_func = -(D16 * wxx+ D26 * wyy + 2* D66 * wxy)
    for x_sub in [0, a]:
            for y_sub in [0, b]:
                M_xy = M_xy_func.subs({x: x_sub, y: y_sub}).evalf()
                M_xy_abs  = max(M_xy_abs, np.abs(M_xy))
                
    
          
    

    
    return deflection, M_x, M_y, M_xy_abs