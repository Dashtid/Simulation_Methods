# Simulation Methods in Biomedical Engineering

This repository contains MATLAB code and scripts for various laboratory assignments in the course **Simulation Methods in Biomedical Engineering**. The course covers a range of numerical methods used in biomedical engineering, including solving partial differential equations (PDEs), boundary value problems, molecular dynamics, and the finite element method (FEM).

## Contents

- **Boundary Value Problems**

  - `Boundary_Value_Problem_Heat_Conduction.m`  
    Solves and visualizes a 1D heat conduction boundary value problem using finite difference methods for different step sizes and convection velocities.

- **Finite Element Method**

  - `1D_FEM_of_2LinearElementSystem.m`  
    Solves a 1D FEM problem using two linear elements and plots the displacement solution.
  - `Finite_Element_Problem_Fishing_Rod.m`  
    Models a fishing rod using the finite element method, assembling global stiffness matrices and plotting the rod's deformation.

- **Molecular Dynamics**
  - `Molecular_Simulation_Verlets_Algorithm.m`  
    Simulates the van der Waals interaction between two Argon atoms using the Lennard-Jones potential and the Velocity Verlet algorithm. Compares results for different time steps.

## How to Use

1. Open the desired `.m` script in MATLAB.
2. Run the script to reproduce the results and plots for each numerical method.
3. Modify parameters within the scripts to experiment with different scenarios.

## Requirements

- MATLAB (R2018b or newer recommended)
- No additional toolboxes are strictly required, but the scripts use standard MATLAB functions for matrix operations and plotting.

## File Overview

- **Boundary_Value_Problem_Heat_Conduction.m**  
  Approximates temperature distribution in a rod for varying mesh sizes and convection velocities.
- **1D_FEM_of_2LinearElementSystem.m**  
  Demonstrates assembling and solving a simple 1D FEM system.
- **Finite_Element_Problem_Fishing_Rod.m**  
  Assembles global stiffness matrices for a multi-segment rod and visualizes its deformation.
- **Molecular_Simulation_Verlets_Algorithm.m**  
  Implements the Velocity Verlet algorithm for molecular dynamics simulation.

## Notes

- Each script is self-contained and includes comments explaining the steps and parameters.
- Example problems are chosen for clarity and educational value; you can adapt the scripts for more complex systems.
- Plots are generated automatically to visualize the results of each simulation.

## License

MIT License

---

_Maintained by David Dashti. For questions or suggestions, please open an issue or contact me._
