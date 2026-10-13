===========================================================================
COMPUTATIONAL FLUID MECHANICS — HOMEWORK 1 
===========================================================================

1. ENVIRONMENT & LANGUAGE VERSION
---------------------------------------------------------------------------
- Language: MATLAB (R2026a)
- Toolboxes / Packages: None required (uses standard MATLAB built-in functions)
- Custom Package: +cfm (included in the root directory)


2. CODE STRUCTURE
---------------------------------------------------------------------------
- run_s1_accuracy.m   : Section 1.1 — Finite-difference accuracy study
- run_s1_convection.m : Section 1.2 — Steady-state 1D convection study
- run_s2_flow.m       : Section 2   — 2D compressible potential flow (FVM)
- +cfm/               : Modular package containing solvers and utilities


3. REPRODUCING FIGURES AND TABLES
---------------------------------------------------------------------------

SECTION 1.1: Accuracy of Finite-Difference Methods
Command:
  run_s1_accuracy
Outputs:
  - Figure 1: Absolute error vs step size h for k = 0.02
  - Figure 2: Absolute error vs step size h for k = 0.2

SECTION 1.2: Steady-State Convection (1D)
Command:
  run_s1_convection
Outputs:
  - figures/s1_convection_solution.pdf    : Analytical vs numerical solution (N=40)
  - figures/s1_convection_error.pdf       : Pointwise local error distribution (N=40)
  - figures/s1_convection_convergence.pdf : L2 and Linf mesh refinement study
  - Command Window                        : Mesh convergence summary table (N=20..320)

SECTION 2: Compressible Potential Flow around a Diamond Airfoil (2D)
Command:
  run_s2_flow

By default (with regime = 'all'), this script runs:
  - Subsonic  (M = 0.3)
  - Transonic (M = 0.91)
  - Supersonic (M = 2.2)

Outputs:
  - figures/s2_flow_M_0.3.pdf             : Subsonic potential, local Mach, density ratio
  - figures/s2_flow_M_0.91.pdf            : Transonic flow field and local supersonic pocket
  - figures/s2_flow_M_2.2.pdf             : Supersonic flow field and Mach wave reflections
  - results/s2_method_boundary_solver.csv : Summary table of iterations, residuals, and solver time

Optional Studies (set in the parameter block at the top of run_s2_flow.m):
  - Set p.DO_VERIF = true:
      Produces: figures/s2_verification_error.pdf
      Prints: Mach-wave angle verification and observed accuracy orders (LaTeX table format)
  - Set p.DO_EXPLO = true:
      Produces: figures/s2_grid_study_cost.pdf
      Prints: Mesh refinement cost vs error table (LaTeX table format)
===========================================================================

**Authorship.** These functions and the corresponding blocks were generated with the help of an AI assistant (Claude, Anthropic, Gemini) and reviewed and run by the group.