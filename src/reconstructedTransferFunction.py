from constants import *
from discreteFilter import *
from filters import *
import matplotlib.pyplot as plt
#-----------------------------------------------------------------------
class ReconstructedTransferFunction:
    """
    Class for defining the reconstructed target transfer function.

    Attributes: 
        discreteFilter (DiscreteFilter): the discrete filter object. 
        coeffs (double []): List holding the discrete filter 
                                coeffs g0, g1, ..., gM. 
        nVanCittert (int): The number of van Cittert iterations
                           (reconstruction order) to calculate
                           the reconstructed target transfer function.

    Author: Z. Nikolaou. 
    """

    def __init__(self, discreteFilter, targetForRecFilter=None, nVanCittert=0):
        self._dfilter=discreteFilter
        if targetForRecFilter is None: 
            self._target=discreteFilter
        else:
            self._target=targetForRecFilter
        if isinstance(self._target, DiscreteFilter):
            isDiscreteTarget=True
        elif isinstance(self._target, Filter):
            isDiscreteTarget=False
        else:
            print('Error: no such filter instance')
            exit()
        self._nV=nVanCittert
        self._forCoeffs=self._dfilter.coeffs
        self._forTfunc=\
                self._dfilter.getTransferFunction(self._forCoeffs)
        if isDiscreteTarget:
            self._recTfunc= \
                self._dfilter.getTransferFunction(self._forCoeffs, \
                self._nV)
        else:
            self._recTfunc=self._target.getTransferFunction(self._nV)

        self.name="Inverse filter at N van Cittert="+str(nVanCittert)

    def name(self):
        return "Inverse filter at N van Cittert="+str(self._nV)

    def getTransferFunction(self): return self._recTfunc

    def getOriginalTransferFunction(self): return self._forTfunc 

    def getVQTransferFunction(self, vCoeffs):
        #returns Vd*Qd
        vgTfunc=self._forTfunc*\
                self._dfilter.getTransferFunction(vCoeffs)
        return vgTfunc

    def getNvanCittert(self): return self._nV
#-----------------------------------------------------------------------
