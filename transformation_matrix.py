import numpy as np

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