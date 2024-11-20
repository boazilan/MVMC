# MVMC

A fully-vectorized MATLAB implementation for Monte Carlo simulations to solve the radiative transfer equation
The statistical approach is known as Monte Carlo for photon transport. The codes provided here are MATLAB implementations corresponding to the C-code based software package MCML, for which the C source codes along with detailed documentation are available on GitHub https://github.com/lhvwang/MCML and references cited there.

These codes use MATLAB's Basic Functions and make use of MATLAB’s vectorized matrix and vector operations. The codes do not loop over the "photons" or use conditional statements. This vectorized approach allows for much faster runtimes using MATLAB, albeit with larger memory requirements.

The codes provided here were used for the computational results in the paper "Asymptotic behavior of the reflectance of a narrow beam by a plane-parallel slab", by Boaz Ilan and Arnold D. Kim, published by the Journal of the Optical Society of Americca A, Vol. 41, No. 12, December 2024.

Directory structure: 
1. Top-level directory: contains most of the main scripts (see below).
2. MVMC: Monte Carlo codes
3. processing: binning and comparison with theory
4. data: results saved .mat files
5. plotting: plotting scripts
6. figures: figures saved as .eps files
   
The main scripts to run are:
1. Compare_MC_First_Order_HG.m: generates Figure 1, which compares  the spatially resolved reflectance with the asymptotic theory using the Henyey-Greenstein phase function.
2. Asymptotic_Error_HG_Vary_Epsilon.m: generates Figure 2, which shows that the asymptotic error scales like O(epsilon^2).
3. Influence_HG_vary_g.m: generates Figure 3, which shows that the influence integral accurately predicts the reflectance using the Henyey-Greenstein phase function.
4. Influence_HG_vary_g.m: generates Figure 4, which shows that the influence integral  accurately predicts the reflectance using the Henyey-Greenstein phase function.
5. plotting/Plot_figure_5.m: generates Figure 5, which shows that the Henyey-Greenstein and TTRM phase functions have completely different behaviors.
6. Reflectance_TTRM_vs_HG_Vary_epsilon.m: generates Figure 6, which compares the reflectances from Monte Carlo simulations using Henyey-Greenstein and TTRM phase functions.
   
Boaz Ilan
