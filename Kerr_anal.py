from einsteinpy.symbolic import (
    MetricTensor,
    ChristoffelSymbols,
    RicciTensor,
    RicciScalar,
    RiemannCurvatureTensor
)
from sympy import symbols, Function, simplify, sin, cos, latex

def compute_and_export_latex(metric, coordinates, filename="tensor_output.tex"):
    dx_ds = symbols('v0 v1 v2 v3')  # Derivatives of coordinates

    metric_tensor = MetricTensor(metric, syms=coordinates)
    christoffel = ChristoffelSymbols.from_metric(metric_tensor)
    ricci_tensor = RicciTensor.from_metric(metric_tensor)
    ricci_scalar = RicciScalar.from_riccitensor(ricci_tensor)
    riemann_tensor = RiemannCurvatureTensor.from_metric(metric_tensor)

    with open(filename, "w") as f:
        f.write("\\documentclass{article}\n")
        f.write("\\usepackage{amsmath}\n\\usepackage{geometry}\n\\geometry{margin=1in}\n")
        f.write("\\begin{document}\n")

        # Christoffel Symbols
        f.write("\\section*{Non-zero Christoffel Symbols}\n\\begin{align*}\n")
        for i in range(metric_tensor.dims):
            for j in range(metric_tensor.dims):
                for k in range(metric_tensor.dims):
                    symbol = simplify(christoffel[i, j, k])
                    if not symbol.is_zero:
                        expr = latex(symbol)
                        f.write(f"\\Gamma^{{{i}}}_{{{j}{k}}} &= {expr} \\\\\n")
        f.write("\\end{align*}\n")

        # Geodesic Equations
        f.write("\\section*{Geodesic Equations}\n\\begin{align*}\n")
        for alpha in range(metric_tensor.dims):
            total = 0
            for beta in range(metric_tensor.dims):
                for gamma in range(metric_tensor.dims):
                    total += christoffel[alpha, beta, gamma] * dx_ds[beta] * dx_ds[gamma]
            simplified = simplify(-total)
            expr = latex(simplified)
            coord_name = latex(coordinates[alpha])
            f.write(f"\\frac{{d^2 {coord_name}}}{{ds^2}} &= {expr} \\\\\n")
        f.write("\\end{align*}\n")

        # Ricci Tensor
        f.write("\\section*{Non-zero Ricci Tensor Components}\n\\begin{align*}\n")
        for i in range(metric_tensor.dims):
            for j in range(metric_tensor.dims):
                rt_component = simplify(ricci_tensor[i, j])
                if not rt_component.is_zero:
                    expr = latex(rt_component)
                    f.write(f"R_{{{i}{j}}} &= {expr} \\\\\n")
        f.write("\\end{align*}\n")

        # Ricci Scalar
        f.write("\\section*{Ricci Scalar}\n")
        scalar_expr = latex(simplify(ricci_scalar.expr))
        f.write(f"\\[ R = {scalar_expr} \\]\n")

        # Riemann Tensor
        f.write("\\section*{Non-zero Riemann Tensor Components}\n\\begin{align*}\n")
        for rho in range(metric_tensor.dims):
            for sigma in range(metric_tensor.dims):
                for mu in range(metric_tensor.dims):
                    for nu in range(metric_tensor.dims):
                        r_component = simplify(riemann_tensor[rho, sigma, mu, nu])
                        if not r_component.is_zero:
                            expr = latex(r_component)
                            f.write(f"R^{{{rho}}}_{{{sigma}{mu}{nu}}} &= {expr} \\\\\n")
        f.write("\\end{align*}\n")

        f.write("\\end{document}\n")

    print(f"LaTeX output written to: {filename}")


if __name__ == "__main__":
    # Coordinates
    t, r, theta, phi = symbols('t r theta phi')
    m = symbols('m')
    c = symbols('c')
    alpha = symbols('alpha')

    
    
    # Spherically symmetric static metric
    g_metric = [
        [(1 - 2*m*r/(r**2 + alpha**2*cos(theta)**2)), B, C, D],
        [E, F, G, H],
        [I, J, K, L],
        [M, N, O, P]
    ]

    coordinates = [t, r, theta, phi]
    compute_and_export_latex(g_metric, coordinates)
