# pyDof

# ![title](figs/graphicFilter.png)
    
**pyDof**, is a library for computing optimised coefficients for discrete forward and inverse filters according to the methodology
developed in [1, 2]. 

The user can select amongst Gaussian, Helmholtz, 3-point Implicit, Butterworth, and may also 
implement their own filter transfer function. The filter coefficients are computed by solving an optimisation problem
and saved to file. Numerous filter constraints are currently supported: positivity, monotonicity etc.    

**Author: Z. Nikolaou (2025).** 

**Contact: ZachariasMNic@gmail.com**

**References:** 

 [1] Z. Nikolaou, L. Vervisch, P. Domingo. Comp. Fluids 255 (2023) 105840. 

 [2] Z. Nikolaou, L. Vervisch, P. Domingo. J. Fluid Mech. 983 (2023) A47. 

Installation: 
-------------

1. Clone the repo to your local machine.  

2. Set environment variables: 
   
   *source setenv.sh*
    
To use:
-------

Have a look at the *./examples/* directory.
