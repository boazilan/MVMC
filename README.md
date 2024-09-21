# MVMC
A fully-vectorized MATLAB implementation for Monte Carlo simulations to solve the radiative transfer equation

The statistical approach is known as Monte Carlo for photon transport. The codes provided here are MATLAB implementations corresponding to the C-code based software package MCML, for which the C source codes along with detailed documentation are available on GitHub https://github.com/lhvwang/MCML and references cited there.


The codes use MATLAB's Basic Functions and make use of MATLAB’s vectorized matrix and vector operations. The codes do not loop over the "photons" or use conditional statements. This vectorized approach allows for much faster runtimes when using MATLAB, albeit with larger memory requirements.

These codes were used for some of results contained in the manuscript "Asymptotic behavior of the reflectance of a narrow beam by a plane-parallel slab" by Boaz Ilan and Arnold D. Kim.

Boaz Ilan
