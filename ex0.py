import numpy as np


def build_vector(v1, v2, v3):
    v = np.array([v1, v2, v3])
    return v


def build_matrix(a,b,c,d,e,f,g,h,i,):
    M = np.array([[a,b,c],
                  [d,e,f], 
                    [g,h,i]])
    return M


def solve_system(n,m, A,B,D):

    #left side n,m 1x6
    left_side = np.concatenate((n, m), axis=0)

    #matrix A,B;B,D 6x6
    matrix = np.block([[A, B], [B, D]])

    
    # Solve the system of equations
    solution = np.linalg.solve(matrix, left_side)
    
    return solution


# # a
# n = build_vector(100, 100, 50)
# m = build_vector(0, 0, 0)

# A = build_matrix(62177, 20391, 0, 20391, 62177, 0, 0,0, 20893)
# B = build_matrix(0, 0, 0, 0, 0, 0, 0, 0,0)
# D = build_matrix(1, 0, 0, 0, 1, 0, 0, 0,1)

# #solving system of equations 6 unknowns

# #left side n,m 1x6
# left_side = np.concatenate((n, m), axis=0)

# #matrix A,B;B,D 6x6
# matrix = np.block([[A, B], [B, D]])

# solution = solve_system(n,m,A,B,D)

# print("Solution a:")
# print(solution)

# #b imax = 3
# imax = 3
# n = build_vector(0, 0, 0)
# m = build_vector(1000, 0, 0)

# D = build_matrix(147.87, 2.81, 0, 2.81, 9.05, 0, 0, 0,3.3)
# D1 = 13/12*(1.27**3)*D
# D = build_matrix(9.05, 2.81, 0, 2.81, 147.87, 0, 0, 0,3.3)
# D2 = 1/12*1.27**3*D
# D3 = D1

# D = D1 + D2 + D3
# for i in range(imax):
#     Ai = build_matrix(1/imax, 0, 0, 0, 1/imax, 0, 0, 0, 1/imax)
#     Bi = build_matrix(0, 0, 0, 0, 0, 0, 0, 0,0)
#     A = A+Ai
#     B = B+Bi
     

# solution = solve_system(n,m,A,B,D)
# print("Solution b:")
# print(solution)

# #c 
# imax = 8

# A = build_matrix(147.870, 2.807, 0, 2.807, 9.053, 0, 0,0, 3.3)
# A1 = 127*A
# A8 = A1

# A = build_matrix(9.053, 2.807, 0, 2.807, 147.870, 0, 0,0, 3.3)
# A2 = 127*A
# A7 = A2

# A = build_matrix(43.9,37.3,34.7,37.3,43.9,34.7,34.7,34.7,37.8)
# A3 = 127*A
# A6 = A3

# A = build_matrix(43.9,37.3,-34.7,37.3,43.9,-34.7,-34.7,-34.7,37.8)
# A4 = 127*A
# A5 = A4
# A = A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8

# B = build_matrix(0, 0, 0, 0, 0, 0, 0, 0,0)

# for i in range(imax):
#     D = build_matrix(1/imax, 0, 0, 0, 1/imax, 0, 0, 0, 1/imax)

# n = build_vector(100, 100, 50)
# m = build_vector(0, 0, 0)

# solution = solve_system(n,m,A,B,D)

# print("Solution c:")
# print(solution)
    
