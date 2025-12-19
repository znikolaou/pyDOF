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
    def aperce(a, bTarget):
        n=len(a)
        nc=0
        e=np.zeros(n)
        for i in range(n):
            if(bTarget[i]!=0.0e0): 
                nc=nc+1
                e[i]=100.0*abs((a[i]-bTarget[i])/bTarget[i])
        return nc, e

    @staticmethod
    def maperc(a, bTarget):
        nc, aperce=Errors.aperce(a, bTarget)
        return np.sum(aperce)/nc

    @staticmethod
    def maxAbsPerc(a, bTarget):
        nc, aperce=Errors.aperce(a, bTarget)
        return [np.argmax(aperce), np.max(aperce)]

    @staticmethod
    def maxOvershoot(a): return (np.max(a)-1.e0)*100.0


#-----------------------------------------------------------------------

