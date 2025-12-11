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

    @staticmethod
    def mperc(a, bTarget):
        n=len(a)
        e=0.0e0
        for i in range(n):
            if(b[i]!=0.0e0): e=e+100.0*(a[i]-b[i])/b[i]
        return e

    @staticmethod
    def maxOvershoot(a): return (np.max(a)-1.e0)*100.0


#-----------------------------------------------------------------------

