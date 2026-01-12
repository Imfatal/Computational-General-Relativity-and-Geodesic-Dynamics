import numpy as np
import sympy as sp
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols
from sympy import symbols, diff, Function

def compute_christoffel_and_geodesics(metric, coordinates):
    # Define metric tensor
    metric_tensor = MetricTensor(metric, syms=coordinates)
    
    # Compute Christoffel symbols
    christoffel = ChristoffelSymbols.from_metric(metric_tensor)

    # Define velocity functions
    tau = symbols("tau")  # Affine parameter (proper time)
    velocity = [Function(f"x{i}")(tau) for i in range(len(coordinates))]

    # Compute geodesic equations
    geodesic_eqs = []
    for i in range(len(coordinates)):
        equation = diff(velocity[i], tau, tau)
        for j in range(len(coordinates)):
            for k in range(len(coordinates)):
                equation += -christoffel[i, j, k] * diff(velocity[j], tau) * diff(velocity[k], tau)
        geodesic_eqs.append(equation)

    return christoffel.tensor(), geodesic_eqs

# Define coordinate symbols (Example: Schwarzschild metric in spherical coordinates)
t, r, theta, phi = symbols("t r theta phi")
coordinates = [t, r, theta, phi]

# Schwarzschild metric (Modify for other metrics)
M = symbols("M")
g_metric = [[-(1 - 2*M/r), 0, 0, 0],
            [0, 1/(1 - 2*M/r), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]]

# Compute Christoffel symbols and geodesic equations
christoffel_symbols, geodesic_equations = compute_christoffel_and_geodesics(g_metric, coordinates)

# Print Christoffel symbols
print("Christoffel Symbols:")
print(christoffel_symbols)

# Print geodesic equations
print("\nGeodesic Equations:")
for eq in geodesic_equations:
    print(eq)
