import math
import numpy as np
#-----------------------------------------------------------------------
class Const:
    """
    Class for holding constants. 

    Author: Z. Nikolaou.
    """
    HALF=0.5e0
    ZERO=0.0e0
    ONE=1.0e0
    TWO=2.0e0
    PI=np.pi
    N=700
    TOL_NYQDAMP=0.1
    TOL_OPTIMISER=1.0e-6
    MAX_ITERATIONS=1500
    #
    NUM_FORMAT_SHORT="{:.8e}.".format
    NUM_FORMAT="{:.16e}".format
    #
    LOGO="PYDOF 24.00"
    LINE="!*************************************************************"
    HEADER_OUT=LINE+"\n! "+LOGO+"\n! Author: Z. Nikolaou"+"\n"+LINE
#-----------------------------------------------------------------------

