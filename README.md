## Blade Element Momentum (BEM) Analysis with Tip Speed Ratio ##

[![MATLAB](https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white)](https://www.mathworks.com/)

Overview
This project implements the Blade Element Momentum (BEM) method to analyze the aerodynamic performance of a horizontal-axis wind turbine.

The code computes:
- Axial induction factor (a)
- Tangential induction factor (a')
- Angle of attack (AoA)
- Flow angle (φ)
- Thrust distribution
- Power output (P)
- Power coefficient (Cp)

Instead of using a fixed tip speed ratio (λ), the code:
- Sweeps λ from **0 to 10**
- Computes Cp at each value
- Finds the **optimal λ (maximum Cp)**
- Plots aerodynamic distributions at that optimal λ condition
---
## Required Files
Ensure the following airfoil data files are present:
```
NACA0018_50k.csv
NACA0018_100k.csv
NACA0018_200k.csv
NACA0018_500k.csv
NACA0018_1000k.csv
```
---
##  Code Workflow

### 1. Input Parameters
- Wind speed (U)
- Blade radius (R)
- Number of blades (B)
- Air density (ρ)
- Viscosity (μ)
- Blade twist and chord distribution

### 2. Blade Discretization
- Blade divided from **20% radius to tip**
- Radial elements used for local calculations

### 3. Airfoil Data Handling
- Uses multiple Reynolds number datasets
- Performs:
  - Interpolation over angle of attack
  - Interpolation between Reynolds numbers

### 4. Iterative BEM Solver
For each radial section:
- Initialize a and a'
- Compute:
  - Flow angle (φ)
  - Angle of attack (α)
  - Lift (Cl) and drag (Cd)
- Apply:
  - **Prandtl tip loss correction**
- Iterate until convergence

### 5. Performance Calculation
- Thrust: dT
- Torque: dQ
- Power: P
- Efficiency: Cp

### 6. Tip Speed Ratio Sweep
- λ varies from **0 → 10**
- Full BEM solution computed at each λ
- Cp stored for all values

### 7. Optimal Condition Selection
- Finds λ with **maximum Cp**
- Stores corresponding aerodynamic results

---
## Outputs
1. Place all CSV files in the same directory as the script  
2. Run the MATLAB script  
3. Observe:
      - Cp vs λ graph  
      - Optimal λ value  
      - Tangential Induction Factor vs Radius
      - AOA Distribution vs Radius
      - Flow Angle vs Radius
      - Thrust Distribution vs Radius
4. Radial Distributions (at optimal λ)
      - Axial induction factor (a)
      - Tangential induction factor (a')
      - Angle of attack (AoA)
      - Flow angle (φ)
      - Thrust distribution (dT)
### Assumptions & Limitations
- No Glauert correction (high induction region)
- No root loss correction
### Future Improvements
- Add Glauert correction
- Include root loss effects

  


