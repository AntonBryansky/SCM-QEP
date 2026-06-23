# Quadratic Eigenvalue Problem formulation of the Spectral Collocation Method
Calculate leaky Lamb wave dispersion and attenuation curves for orthotropic layers bonded to isotropic substrates along principal axes of symmetry.

## Description
The MATLAB-based script for calculation of dispersion and attenuation curves for orthotropic plates bonded to isotropic substrate in decoupled Lamb leaky wave mode using Spectral Collocation Method.

## Paper
A. A. Bryansky, R. R. Badykov, Paper name, Journal name (year) Volume Pages-Pages.

## Files

`SCM_SL.m` and `SCM_Lamb_SL.m` - Spectral Collocation Method for calculation of the dispersion curves for free orthotropic single-layer media, coupled and Lamb decouplde cases, respectively.

`SCM_ML.m` and `SCM_Lamb_ML.m` - Spectral Collocation Method for calculation of the dispersion curves for free orthotropic multi-layer media, coupled and Lamb decouplde cases, respectively.

`SCM_QUADEIG_SL.m` and `SCM_Lamb_QUADEIG_SL.m` - Spectral Collocation Method reformulated as Quadratic Eigenvalue Problem for calculation of the dispersion curves for free orthotropic single-layer media, coupled and Lamb decouplde cases, respectively.

`SCM_QUADEIG_ML.m` and `SCM_Lamb_QUADEIG_ML.m` - Spectral Collocation Method reformulated as Quadratic Eigenvalue Problem for calculation of the dispersion curves for free orthotropic multi-layer media, , coupled and Lamb decouplde cases, respectively.

`SCM_Lamb_QUADEIG_SL_SUB_POTENTIAL.m` - Spectral Collocation Method reformulated as Quadratic Eigenvalue Problem for calculation of the dispersion curves for the orthotropic single-layered media bonded to the isotropic substrate.

`balance2.m` and `gcg.m` - functions required for balancing algorithm (look [Important links](#important-links)).

`quadeig.m` - Quadratic eigenvalue problem solver (look [Important links](#important-links))

## Literature
1. F. Hernando Quintanilla, M. J. S. Lowe, R. V. Craster, Modeling guided elastic waves in generally anisotropic media using a spectral collocation method, J. Acoust. Soc. Am. 137.3 (2015) 1180-1194. [DOI](https://doi.org/10.1121/1.4913777)
2. F. Hernando Quintanilla, M. J. S. Lowe, R. V. Craster, The symmetry and coupling properties of solutions in general anisotropic multilayer waveguides, J. Acoust. Soc. Am. 141(1) (2017) 406-418. [DOI](https://doi.org/10.1121/1.4973543)
3. M. Mekkaoui, S. Nissabouri, H. Rhimini, Towards an Optimization of the Spectral Collocation Method with a New Balancing Algorithm for Plotting Dispersion Curves of Composites with Large Numbers of Layers, J. Appl. Comput. Mech. 10(4) (2024) 801-816. [DOI](https://doi.org/10.22055/jacm.2024.45578.4390)
4. I. Zitouni, H. Rhimini, A. Chouaf, Comparative study of the spectral method, DISPERSE and other‎ classical methods for plotting the dispersion curves in‎ anisotropic plates, J. Appl. Comput. Mech. 9(4) (2023) 955-973. [DOI](https://doi.org/10.22055/jacm.2023.42530.3941)
5. Georgiades, E., Lowe, M. J., & Craster, R. V. Leaky wave characterisation using spectral methods, J. Acoust. Soc. Am. 152(3) (2022) 1487-1497. [DOI](https://doi.org/10.1121/10.0013897)
6. Georgiades, E., Lowe, M. J., & Craster, R. V. Computing leaky Lamb waves for waveguides between elastic half-spaces using spectral collocation, J. Acoust. Soc. Am. 155(1) (2024) 629-639. [DOI](https://doi.org/10.1121/10.0024467)
7. Hammarling, S., Munro, C. J., & Tisseur, F. An algorithm for the complete solution of quadratic eigenvalue problems. ACM Transactions on Mathematical Software (TOMS) 39(3) (2013) 1-19. [DOI](https://doi.org/10.1145/2450153.2450156)
8. A. Huber, M. G. Sause, Classification of solutions for guided waves in anisotropic composites with large numbers of layers, J. Acoust. Soc. Am. 144(6) (2018) 3236-3251. [DOI](https://doi.org/10.1121/1.5082299)
9. J. Weideman, S. Reddy, A MATLAB differentiation matrix suite, ACM Trans. Math. Softw. 26(4) (2000) 465–519. [DOI](https://doi.org/10.1145/365723.365727)

## Important links:
[Dispersion Calculator](https://github.com/ArminHuber/Dispersion-Calculator) by Armin Huber

About algorithm for solving of quadratic eigenvalue problem:
[quadeig](https://github.com/ftisseur/quadratic-eigensolver) 

About balancing algorithm for eigenvalue problem:
[BALANCE2 Balancing generalized eigenvalue problem](https://www.mathworks.com/matlabcentral/fileexchange/49719-balance2-balancing-generalized-eigenvalue-problem) (requires
[GCG Generalized conjugate gradient method](https://www.mathworks.com/matlabcentral/fileexchange/49720-gcg-generalized-conjugate-gradient-method))

## Thanks to
* Dr. Armin Huber, German Aerospace Center (DLR), Augsburg, Germany
* [asm-jaime](https://github.com/asm-jaime)
