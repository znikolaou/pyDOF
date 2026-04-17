from utils import *
from filters import *
from discreteFilter import *
from errors import *
from reconstructedTransferFunction import *
#-----------------------------------------------------------------------
class ObjectiveFunction(Errors, Const):
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
        self._trgFunc=self._target.getTransferFunction() 
        self._counter=0
        self._lastEr=None
        self._lastSol=None
        self._lamda=lamda
        self._er=error
        
    def get(self, coeffs):
        self._counter+=1
        gradEr=Const.ZERO
        grad=Const.ZERO

        if(isinstance(self._target, ReconstructedTransferFunction)):
            self._lastSol=self._target.getVQTransferFunction(coeffs)
        else:
            self._lastSol=self._dfilter.getTransferFunction(coeffs)
        
        if(self._lamda!=Const.ZERO): grad=np.gradient(self._lastSol)
 
        if self._er=='msqe':
            self._lastEr=self.msqe(self._trgFunc, self._lastSol)+self._lamda*self.meanSq(grad)
        elif self._er=='mabse':
            self._lastEr=self.mabse(self._trgFunc, self._lastSol)+self._lamda*self.meanAbs(grad)
 
        return self._lastEr

    def getLastSolution(self): return self._lastSol

    def getLastError(self): return self._lastEr

#-----------------------------------------------------------------------

