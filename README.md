# Quasi-1D CD Nozzle — Steger-Warming FVS Solver

Numerical solver for compressible, inviscid, quasi-1D flow through a convergent-divergent nozzle, implemented as the final project of the **Computational Fluid Dynamics** course (ESZS 035-17 — 2025.2), Universidade Federal do ABC.

The solver uses the **Steger-Warming Flux Vector Splitting** scheme on a finite volume framework and is validated against the isentropic analytical solution, with automatic detection of the flow regime (subsonic, supersonic, or internal normal shock).

---

## Problem Description

The governing equations are the quasi-1D Euler equations in conservative form:

$$\frac{\partial \mathbf{U}}{\partial t} + \frac{1}{S}\frac{\partial (S\mathbf{F})}{\partial x} = \mathbf{Q}$$

where the state, flux, and source term vectors are:

$$\mathbf{U} = \begin{bmatrix} \rho \\ \rho u \\ e \end{bmatrix}, \quad \mathbf{F} = \begin{bmatrix} \rho u \\ \rho u^2 + p \\ (e+p)u \end{bmatrix}, \quad \mathbf{Q} = \begin{bmatrix} 0 \\ \dfrac{p}{S}\dfrac{dS}{dx} \\ 0 \end{bmatrix}$$

The system is closed by the equation of state for a calorically perfect gas:

$$p = \left(e - \frac{\rho u^2}{2}\right)(\gamma - 1)$$

### Nozzle Geometry

The cross-sectional area varies parabolically along the nozzle axis:

$$S(x) = 0.8x^2 - 0.6x + 0.35, \quad x \in [0,\, 1] \text{ m}$$

The throat is located at $x^* = 0.375$ m, where $S^* = 0.2375$ m².

### Non-Dimensionalization

All variables are non-dimensionalized using the following reference state:

| Quantity | Value |
|---|---|
| $\rho_\text{ref}$ | 1.486 kg/m³ |
| $a_\text{ref}$ | 320.44 m/s |
| $l_\text{ref}$ | 1 m |
| $T_\text{ref}$ | 357.778 K |
| $R_\text{ref}$ | 287 J/(kg·K) |

---

## Numerical Method

### Finite Volume Discretization

The domain is discretized using a **cell-vertex** scheme with **400 equally spaced nodes**. Half-volumes are used at the inlet and outlet boundaries. The explicit Euler time-marching scheme gives:

$$\mathbf{U}_i^{(n+1)} = \mathbf{U}_i^{(n)} - \frac{\Delta t}{\Delta x \, \bar{S}_i}\left[S_{i+1/2}\mathbf{F}_{i+1/2}^{(n)} - S_{i-1/2}\mathbf{F}_{i-1/2}^{(n)}\right] + \mathbf{Q}_i^{(n)}\Delta t$$

where $S_{i\pm 1/2} = S\!\left(\frac{x_i + x_{i\pm 1}}{2}\right)$ and $\bar{S}_i = \frac{S_{i+1/2} + S_{i-1/2}}{2}$.

### Steger-Warming Flux Vector Splitting

The flux Jacobian $\mathbf{A} = \partial\mathbf{F}/\partial\mathbf{U}$ is decomposed as $\mathbf{A} = \mathbf{A}^+ + \mathbf{A}^-$, where:

$$\mathbf{A}^\pm = \mathbf{P}\mathbf{\Lambda}^\pm\mathbf{P}^{-1}, \qquad \lambda^\pm = \frac{\lambda \pm |\lambda|}{2}$$

with eigenvalues $\lambda \in \{u,\; u+a,\; u-a\}$. The interface fluxes are:

$$\mathbf{F}_{i+1/2} = \mathbf{A}^+_{(i)}\mathbf{U}_i + \mathbf{A}^-_{(i+1)}\mathbf{U}_{i+1}$$
$$\mathbf{F}_{i-1/2} = \mathbf{A}^-_{(i)}\mathbf{U}_i + \mathbf{A}^+_{(i-1)}\mathbf{U}_{i-1}$$

### CFL Condition

The time step is controlled via:

$$\Delta t = \text{CFL} \cdot \frac{\Delta x}{\max_i |u_i + a_i|}, \qquad \text{CFL} = 0.5$$

---

## Boundary Conditions

Boundary conditions are imposed via Riemann invariants (method of characteristics).

### Inlet — Subsonic

Two characteristics enter the domain. The stagnation pressure $p_0$ and stagnation temperature $T_0$ are fixed. The free degree of freedom is resolved through the outgoing characteristic (energy equation), yielding $\delta u_1$, from which $p_1$, $T_1$, $\rho_1$, and $e_1$ are updated via isentropic relations.

### Outlet — Subsonic

Static pressure $p_\text{out}$ is fixed ($\delta p = 0$). Density and velocity are updated from the two outgoing characteristics (mass and momentum equations).

### Outlet — Supersonic

All characteristics exit the domain. Velocity, pressure, and density are all extrapolated from the interior using the three characteristic equations.

---

## Simulated Cases

Both cases use air as the working fluid ($\gamma = 1.4$, $R = 287$ J/(kg·K)) with identical inlet initial conditions.

|  | **Case 1** | **Case 2** |
|---|---|---|
| $M_\text{in}$ | 0.7 | 0.7 |
| $p_\text{in}$ | 108.99 kPa | 108.99 kPa |
| $T_\text{in}$ | 255.56 K | 255.56 K |
| $p_\text{out}$ | **85.00 kPa** | **40.00 kPa** |
| Expected regime | Subsonic / mild | Normal shock inside nozzle |

---

## Automatic Regime Detection

The post-processing script identifies the flow regime from the converged numerical solution before computing the analytical reference:

| Condition | Detected Regime |
|---|---|
| $M_\text{exit} > 1$ | Fully supersonic |
| $M_\text{max,div} > 1.05$ **and** $\Delta p_0/p_{0,\text{in}} > 0.5\%$ | Normal shock inside nozzle |
| Otherwise | Fully subsonic |

For the **internal shock** case, the shock position is estimated from the steepest negative $\nabla M$ in the divergent section. The post-shock analytical solution is computed using Rankine-Hugoniot relations with the corrected stagnation pressure $p_{0,2} = p_{0,1} \cdot (p_{02}/p_{01})$.

---

## Repository Structure

```
.
├── steger_warming.cpp   # Main solver: time marching, BCs, file output
├── functions.cpp        # Helpers: linspace, nozzle area, FVS matrices, residuals
├── functions.h          # Header file
├── plot.py              # Post-processing: comparison plots and error analysis
└── README.md
```

### Output files (generated by the solver)

| File | Contents |
|---|---|
| `pressao.txt` | Static pressure $p(x)$ |
| `pressao_estag.txt` | Stagnation pressure $p_0(x)$ |
| `temperatura.txt` | Static temperature $T(x)$ |
| `temp_estag.txt` | Stagnation temperature $T_0(x)$ |
| `mach.txt` | Mach number $M(x)$ |
| `residuo1.txt` | Density residual history |
| `residuo2.txt` | Momentum residual history |
| `residuo3.txt` | Energy residual history |

---

## How to Run

### 1. Compile and run the solver

```bash
g++ -O2 -std=c++17 -o solver steger_warming.cpp functions.cpp
./solver
```

### 2. Run the post-processing script

```bash
pip install numpy matplotlib scipy
python plot.py
```

The script generates:

| File | Contents |
|---|---|
| `comparacao_analitica.png` | Numerical vs. analytical — 2×2 panel ($M$, $p$, $T$, $p_0$) |
| `erro_analitico.png` | Pointwise relative error — 2×2 panel |
| `mach.png` | Individual Mach number plot |
| `pressao.png` | Individual static pressure plot |
| `pressao_estag.png` | Individual stagnation pressure plot |
| `temperatura.png` | Individual static temperature plot |
| `temp_estag.png` | Individual stagnation temperature plot |
| `residuos.png` | Convergence history of all three residuals |

---

## Validation

The analytical solution is computed by solving the isentropic area-Mach relation via bisection (`scipy.optimize.brentq`) on the subsonic and supersonic branches. In isentropic regions, the relative error between the numerical and analytical solutions is verified for $M$, $p$, $T$, and $p_0$ using both $L^\infty$ (maximum) and $L^2$ norms. The stagnation pressure $p_0$ is used as a consistency indicator: it must remain constant throughout the domain in any isentropic region.

---

## Dependencies

| Tool | Minimum version |
|---|---|
| g++ | 11 (C++17) |
| Python | 3.9 |
| NumPy | 1.23 |
| Matplotlib | 3.6 |
| SciPy | 1.9 |

---

## References

- ANDERSON, J. D. *Modern Compressible Flow: With Historical Perspective*. 3rd ed. McGraw-Hill, 2003.
- HIRSCH, C. *Numerical Computation of Internal and External Flows, Vol. 1*. John Wiley & Sons, 1988.
- HIRSCH, C. *Numerical Computation of Internal and External Flows, Vol. 2*. John Wiley & Sons, 1990.
- STEGER, J. L.; WARMING, R. F. Flux vector splitting of the inviscid gasdynamic equations with applications to finite difference methods. *Journal of Computational Physics*, v. 40, p. 263–293, 1981.
