from utils import *
from filters import *
from discreteFilter import *
from errors import *
from reconstructedTransferFunction import *
import matplotlib.pyplot as plt 
#-----------------------------------------------------------------------
class ObjectiveFunction(Errors):
    """
    Class for calculating the optimiser objective function. 

    Attributes: 
        target (Filter or ReconstructedTransferFunction): Target. 
        discreteFilter (DiscreteFilter): 
            The discrete filter object whose coeffs must be found.
        lamda (double): regularisation/smoothness parameter. 

    Author: Z. Nikolaou.
    """ 

    def __init__(self, target, discreteFilter, lamda=0.0e0, error='msqe'):
        self._target=target
        self._dfilter=discreteFilter
        self._k=self._dfilter.getK()
        self._targetTfunc=self._target.getTransferFunction() 
        self._counter=0
        self._lastError=None
        self._lastSolution=None
        self._lamda=lamda
        self._error=error
        
    def get(self, coeffs):
        self._counter+=1
        if(isinstance(self._target, ReconstructedTransferFunction)):
            self._lastSolution=self._target.getVQTransferFunction(coeffs)
        else:
            self._lastSolution=self._dfilter.getTransferFunction(coeffs)
        
        gradEr=0.0e0
        if(self._lamda!=0.0e0): 
            grad=np.gradient(self._lastSolution)
            if self._error=='msqe': 
                gradEr=self._lamda*self.meanSq(grad)
            elif self._error=='mabse':
                gradEr=self._lamda*self.meanAbs(grad)
 
        if self._error=='msqe':
            er=self.msqe(self._targetTfunc, self._lastSolution)
        elif self._error=='mabse':
            er=self.mabse(self._targetTfunc, self._lastSolution)

        self._lastError=er+gradEr
        
        return self._lastError

    def getLastSolution(self): return self._lastSolution

    def getLastError(self): return self._lastError

#-----------------------------------------------------------------------

