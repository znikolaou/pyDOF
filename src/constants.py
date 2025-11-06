import numpy as np
import datetime
#-----------------------------------------------------------------------
class Const:
    """
    Class for holding constants. 

    Author: Z. Nikolaou.
    """

    ZERO=0.0e0
    HALF=0.5e0
    ONE=1.0e0
    TWO=2.0e0
    PI=np.pi
    NK_NOMINAL=400
    DIVERGENCE_FACTOR=1.10
    TOL_NYQDAMP=0.1
    TOL_OPTIMISER=1.0e-6
    MAX_ITERATIONS=500
    #
    NUM_FORMAT_SHORT="{:.8e}.".format
    NUM_FORMAT_LONG="{:.16e}".format
    #
    DATE=datetime.datetime.now()
    MONTH=DATE.month
    YEAR=DATE.year
    #
    LOGO="! PYDOF 2025"
    AUTHOR='! Author: Z. Nikolaou '
    IS_AUTO='! This is an auto-generated file: '+str(DATE)
    LINE="!-------------------------------------------------------------"
    HEADER_OUT=LINE+"\n"+LOGO+"\n"+IS_AUTO+"\n"+AUTHOR+"\n"+LINE
    #-----------------------------------------------------------------------

