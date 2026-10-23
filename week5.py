#Exercise 4: Vary the load with x and find the load factor at which each failure criterion first predicts failure (middle of the ply) with 5 different failure criteria
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

# Find the load factor x at which each failure criterion first predicts failure.
# x = 1 corresponds to the reference loads defined below.
N_per_x = np.array([300.0, -50.0, 25.0])  # forces in N/mm all loads are scaled with x
M_per_x = np.array([50.0, 0.0, 0.0])  # moments in N all loads are scaled with x

###FUNCTIONS##################################################################################

def local_stresses_at_load(load_factor):
    """Return local stresses at the bottom, middle, and top of every ply dependent on load_factor x."""
    strains, curvatures = solve_system(
        load_factor * N_per_x,
        load_factor * M_per_x,
        A,
        B,
        D,
    )

    stresses = np.zeros((len(laminate) * 3, 3))
    for i, theta in enumerate(laminate):
        z_coords = [z[i], (z[i] + z[i + 1]) / 2, z[i + 1]]
        for j, z_coord in enumerate(z_coords):
            stresses[i * 3 + j] = find_local_stress(
                strains, curvatures, z_coord, theta, Q_local
            )
    return stresses


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


def first_failure_load(criterion, location=3, ply_number=None, tolerance=1e-8):
    """Find the first failure at a selected ply location.

    Set location to 0 for bottom, 1 for middle, 2 for top, or 3 for any.
    """
    if location not in (0, 1, 2, 3):
        raise ValueError("location must be 0 (bottom), 1 (middle), 2 (top), or 3 (any)")
    if ply_number is not None and not 1 <= ply_number <= len(laminate):
        raise ValueError("ply_number must be between 1 and the number of plies")

    def has_failed(load_factor):
        stresses = local_stresses_at_load(load_factor)
        failure_flags = criterion(stresses).reshape(len(laminate), 3)
        if ply_number is not None:
            failure_flags = failure_flags[ply_number - 1:ply_number]
        if location == 3:
            return np.any(failure_flags == 1)
        return np.any(failure_flags[:, location] == 1)

    lower = 0.0
    upper = 1.0
    while not has_failed(upper):
        upper *= 2.0

    while upper - lower > tolerance * max(1.0, upper):
        midpoint = (lower + upper) / 2.0
        if has_failed(midpoint):
            upper = midpoint
        else:
            lower = midpoint

    stresses = local_stresses_at_load(upper)
    failure_flags = criterion(stresses).reshape(len(laminate), 3)
    if ply_number is not None:
        selected_flags = failure_flags[ply_number - 1:ply_number]
    else:
        selected_flags = failure_flags

    if location == 3:
        failed_points = np.argwhere(selected_flags == 1)
        selected_ply = int(failed_points[0][0])
        selected_location = int(failed_points[0][1])
    else:
        failed_plies = np.flatnonzero(selected_flags[:, location] == 1)
        selected_ply = int(failed_plies[0])
        selected_location = location

    if ply_number is None:
        ply_number = selected_ply + 1
    point_name = ("bottom", "middle", "top")[selected_location]
    return upper, ply_number, point_name
# Build the Q_matrix
Q_local = ply(E1, E2, nu12, G12)

# Build the A,B,D matrices and coordinate system for the laminates
A, B, D, z= laminate_properties(laminate, t, Q_local)

# Use 0 for bottom, 1 for middle, 2 for top, or 3 to check all locations.
location_to_check = 1
failure_loads = []
for criterion_name, criterion in criteria.items():
    for ply_number in range(1, len(laminate) + 1):
        load_factor, _, point_name = first_failure_load(
            criterion,
            location=location_to_check,
            ply_number=ply_number,
        )
        failure_loads.append(
            {
                "criterion": criterion_name,
                "ply": f"ply {ply_number}",
                "failure_load_factor_x": load_factor,
                "location": point_name,
            }
        )

failure_loads_df = pd.DataFrame(failure_loads)
failure_factor_table = failure_loads_df.pivot(
    index="criterion", columns="ply", values="failure_load_factor_x"
)
print("\nFirst failure load for every ply:")
print(failure_factor_table.to_string(float_format=lambda value: f"{value:.3g}"))

