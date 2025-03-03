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
    N=200
    TOL_NYQDAMP=0.1
    TOL_OPTIMISER=1.0e-6
    MAX_ITERATIONS=500
    #
    NUM_FORMAT_SHORT="{:.8e}.".format
    NUM_FORMAT="{:.16e}".format
    #
    DATE=datetime.datetime.now()
    MONTH=DATE.month
    YEAR=DATE.year
    #
    LOGO="PYDOF "+str(YEAR)+'.'+str(MONTH).zfill(2)
    LINE="!*************************************************************"
    HEADER_OUT=LINE+"\n! "+LOGO+"\n! Author: Z. Nikolaou"+"\n"+LINE
    #
    FORTRAN_HEADER='pydof_common.h'
    F90MODULE='pydofFilterDataModule.f90'
#-----------------------------------------------------------------------

