## Quasi-1D Compressible Flow through Converging-Diverging Rocket Nozzles

This repository contains the numerical implementation and technical report for analyzing quasi-one-dimensional, unsteady, compressible, inviscid fluid flow through a Converging-Diverging (CD) rocket nozzle. The governing time-dependent Euler equations in conservative form are discretized using MacCormack's explicit finite difference scheme.

---

## Overview

The solver simulates steady-state solutions across two primary flow regimes by marching in time:
1. **Isentropic Subsonic Flow**
2. **Isentropic Subsonic-Supersonic (Choked) Flow**

Numerical results for pressure ratio, density ratio, temperature ratio, and Mach number profiles are validated against both analytical compressible flow equations and NASA/NPARC experimental validation data.

---

## Key Features

* Custom-coded MacCormack Scheme (Predictor-Corrector steps with 2nd-order spatial and temporal accuracy).
* Native implementation of algebraic system solvers and array manipulations without relying on high-level CFD library functions.
* Characteristic-based boundary condition enforcement based on local eigenvalue signatures ($\lambda = V, V \pm a$).
* Dynamic CFL-based stability calculation for time-step determination ($\Delta t$).
* Comparative post-processing against analytical functions and NASA CDV experimental benchmark data.


