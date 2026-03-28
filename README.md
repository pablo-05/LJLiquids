# Calculation of g(r) for a Lennard-Jones Fluid

## Project Description
This project is part of an Advanced Statistical Physics course. Its main objective is to calculate the radial distribution function, g(r), for a Lennard-Jones fluid. This is achieved through two distinct approaches: a theoretical resolution using iterative methods and a 3D computational simulation.

## Theoretical Framework
The theoretical section solves the liquid state equations using the following elements:
* **Interaction Potential:** The Lennard-Jones potential is used. Distance r is taken in units of sigma, and epsilon is set to 1. For small values of r, the exponential of the potential is manually set to zero to avoid divergences.
* **Fundamental Equation:** The starting point is the Ornstein-Zernike (OZ) equation.
* **Closure Relations:** Both Percus-Yevick (PY) and Hypernetted-Chain (HNC) closures are implemented.
* **Iterative PY Resolution:** A function f(r) = h(r) - c(r) is defined and iterated using a mixing parameter (alpha) to ensure convergence.
* **Iterative HNC Resolution:** The HNC closure modifies the process to iterate directly using the function h(r).

## Monte Carlo Simulation
The computational simulation adapts standard 2D Monte Carlo techniques for hard spheres and Lennard-Jones fluids into 3D.
* **Hard Spheres:** The base model places hard spheres in a cubic cell with periodic boundary conditions. The maximum random displacement is dynamically adjusted to keep the movement acceptance rate around 0.5.
* **Lennard-Jones:** The final version incorporates the Lennard-Jones potential. Proposed movements are accepted or rejected by evaluating the new energy of the system.
* **Measurement of g(r):** Measurements of the radial distribution function are performed at regular intervals once the system reaches thermal equilibrium.

## Numerical Methods and Transforms
* The theoretical calculations rely heavily on 3D radial Fourier transforms.
* The numerical implementation uses standard algorithms (such as sinft, four1, and realft) for these transformations.

## Execution Considerations
* Care must be taken to handle potential division-by-zero errors at r=0 during 3D inverse transforms (e.g., using the second index instead of the first).
* To achieve convergence at high densities, calculations should start at low densities and increase progressively.
