# pyDof

![title](figs/2dimages_c.png)
![title](figs/2dimages_b.png)
![title](figs/2dimages_e.png)
![title](figs/2dimages_turbvflame.png) 
        
**pyDOF** is a library for designing discrete/digital symmetric filters using 
a multi-objective and constrained optimisation framework [1, 2]. 
This framework allows the user to impose a wide range of constraints to the
discrete filter transfer function such as  monotonicity, positivity, value-fixing, 
gradient smoothing etc. The filter stencil size can be set by the user or computed dynamically
by the solver. **pyDOF** solves the constrained optimisation problem, and outputs automatically 
the discrete filter coefficients in a plain text file which can be readily parsed by an programming language. 
The user can select amongst a wide range of built-in target trasfer functions 
(Gaussian, Helmholtz, 3-point Implicit, Butterworth, sharp cut-off), and may also 
implement their own target transfer function. **pyDOF** can be used to design a wide range
of filters such as low-pass, high-pass, multi band-pass/band-stop etc. 

**Author: Z. Nikolaou (2025).** 

**Contact: ZachariasMNic@gmail.com**

Installation: 
-------------

1. Clone the repo to your local machine.  

   E.g. with ssh: 

   git clone git@github.com:znikolaou/pyDOF.git

2. Navigate to the cloned directory and install required libraries:

   pip install -r requirements.txt

3. Set the environment variables: 
   
   *source setenv.sh*

   This will set up the following environment variables: 

   PYDOF\_HOME : installation directory

   PYDOF\_SRC : source file directory

   PYDOF\_EXAMPLES : examples directory

   PYDOF\_TEST : testing directory

To use:
-------

Navigate to your preferred working directory, and 
add the PYDOF\_SRC environment variable to your Python PATH
so that all the classes under PYDOF\_SRC are available. 
We suggest creating an include.py file and adding: 

import sys
sys.path.append(os.environ.get("PYDOF_SRC"))
from plotlib import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import *
from errors import *

Then in your main.py file simply import the above:

from include import *

Have a look at the *./examples/* directory where 
a sample include.py file is given-this file is then
imported for every example therein as detailed above.

Cite:
-----

 [1] Z. Nikolaou, P. Domingo, L. Vervisch, D. Drikakis. [arXiv link](https://arxiv.org/abs/2606.26830)

 pyDOF: a Python library for the design of discrete forward and inverse filters.

 [2] Z. Nikolaou, L. Vervisch, P. Domingo. Comp. Fluids 255 (2023) 105840. 

 An optimisation framework for the development of explicit discrete forward and inverse filters.

 [3] Z. Nikolaou, L. Vervisch, P. Domingo. J. Fluid Mech. 983 (2023) A47. 

 Revisiting the modelling framework for the unresolved scalar variance.
