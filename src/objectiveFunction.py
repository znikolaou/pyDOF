from utils import *
from filters import *
from discreteFilter import *
from errors import *
from reconstructedTransferFunction import *
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

    def __init__(self, target, discreteFilter, lamda=0.0e0):
        self._target=target
        self._dfilter=discreteFilter
        self._targetTfunc=self._target.getTransferFunction() 
        self._counter=0
        self._lastError=None
        self._lastSolution=None
        self._lamda=lamda
        
    def get(self, coeffs):
        self._counter+=1
        if(isinstance(self._target, ReconstructedTransferFunction)):
            self._lastSolution=self._target.getVQTransferFunction(coeffs)
        else:
            self._lastSolution= self._dfilter.getTransferFunction(coeffs)
           
        self._lastError=self.msqe(self._targetTfunc, self._lastSolution)
        self._gradError=\
                np.mean(np.abs(np.gradient(self._targetTfunc)-np.gradient(self._lastSolution)))

        self._lastError=self._lastError+self._lamda*self._gradError

        return self._lastError

    def getLastSolution(self): return self._lastSolution

    def getLastError(self): return self._lastError

#-----------------------------------------------------------------------

