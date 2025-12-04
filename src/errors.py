import numpy as np
import scipy
#-----------------------------------------------------------------------
class Errors:
    """
    Class for defining error measures for optimiser. 

    Author: Z. Nikolaou
    """

    @staticmethod
    def msqe(a, b): return np.square(np.subtract(a, b)).mean()

    @staticmethod
    def mase(a, b): return np.abs(np.subtract(a, b)).mean()

    @staticmethod
    def msqeSimpson(x, a, b):
        l=np.max(x)-np.min(x)
        return 1.0/l*scipy.integrate.simpson(np.square(np.subtract(a, b)), x)

    @staticmethod
    def msqeTrapz(x, a, b):
        l=np.max(x)-np.min(x)
        return 1.0/l*np.trapz(np.square(np.subtract(a, b)), x)


#-----------------------------------------------------------------------

