from constants import *
from discreteFilter import *
from filters import *
from utils import *
from scipy.signal import savgol_filter
#-----------------------------------------------------------------------
class ReconstructedTransferFunction(Log):
    """
    Class for defining the reconstructed target transfer function.

    Attributes: 
        discreteFilter (DiscreteFilter): the discrete filter object. 
        targetForRecFilter (DiscreteFilter or Filter): target for the
            reconstructed transfer function.
        coeffs (double []): List holding the discrete filter 
                            coeffs g0, g1, ..., gM. 
        nVanCittert (int): The number of van Cittert iterations
                           (reconstruction order) to calculate
                           the reconstructed target transfer function.

    Author: Z. Nikolaou. 
    """

    def __init__(self, discreteFilter, targetForRecFilter=None, nVanCittert=0):
        Log.__init__(self, type(self).__name__)
        self._dfilter=discreteFilter
        self._target=None
        self.__setTarget(targetForRecFilter)
        self._isDiscreteTarget=None
        self.__setIsDiscreteTarget()
        self._nV=nVanCittert
        self._forCoeffs=self._dfilter.coeffs
        self._forTfunc=self._dfilter.getTransferFunction(self._forCoeffs)
        self._recTfunc=None
        self.__setRecTfunc()
        self.name="Inverse filter at N van Cittert="+str(nVanCittert)
    
    def __setTarget(self, targetForRecFilter):
        if targetForRecFilter is None: 
            self._target=self._dfilter
        else: 
            self._target=targetForRecFilter
        return

    def __setIsDiscreteTarget(self):
        if isinstance(self._target, DiscreteFilter): 
            self._isDiscreteTarget=True
        elif isinstance(self._target, Filter): 
            self._isDiscreteTarget=False
        else:
            self.printError('__setIsDiscreteTarget', 'No such filter instance')
            exit()
        return

    def __setRecTfunc(self):
        if self._isDiscreteTarget:
            self._recTfunc=self._dfilter.getTransferFunction(self._forCoeffs, self._nV)
        else:
            self._recTfunc=self._target.getTransferFunction(self._nV)
        return

    def name(self): return "Inverse filter at N van Cittert="+str(self._nV)

    def getTransferFunction(self): return self._recTfunc

    def getOriginalTransferFunction(self): return self._forTfunc 

    def getVQTransferFunction(self, vCoeffs):
        return self._forTfunc*self._dfilter.getTransferFunction(vCoeffs)

    def getNvanCittert(self): return self._nV
#-----------------------------------------------------------------------
