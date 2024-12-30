from utils import *
from filters import *
from discreteFilter import *
from errors import *
import matplotlib.pyplot as plt
from reconstructedTransferFunction import *
#-----------------------------------------------------------------------
class ObjectiveFunction(Errors):
    """
    Class for the optimiser objective function. 

    Attributes: 
        targetOb (Filter, ReconstructedTransferFunction): Target. 
        discreteFilterOb (DiscreteFilter): The discrete filter object.

    Author: Z. Nikolaou.
    """ 
    def __init__(self, target, discreteFilter):
        self._target=target
        self._dfilter=discreteFilter
        self._targetTfunc=self._target.getTransferFunction() 
        self._counter=1
        self._lastError=None
        self._lastSolution=None
        
    def get(self, coeffs):
        self._counter+=1
        if(isinstance(self._target, ReconstructedTransferFunction)):
            self._lastSolution= \
                    self._target.getVQTransferFunction(coeffs)
        else:
            self._lastSolution= \
                    self._dfilter.getTransferFunction(coeffs)
        self._lastError=self.msqe(self._targetTfunc, self._lastSolution)
        return self._lastError

    def getLastSolution(self): return self._lastSolution

    def getLastError(self): return self._lastError
#-----------------------------------------------------------------------

