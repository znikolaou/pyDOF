from constants import *
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

    def __init__(self, discreteFilter, nVanCittert):
        self._dfilter=discreteFilter
        self._forCoeffs=self._dfilter.coeffs
        self._nV=nVanCittert
        self._forTfunc=\
                self._dfilter.getTransferFunction(self._forCoeffs)
        self._recTfunc= \
                self._dfilter.getTransferFunction(self._forCoeffs, \
                nVanCittert)
        self.name="Inverse filter at N van Cittert="+str(nVanCittert)

    def name(self):
        return "Inverse filter at N van Cittert=" \
                +str(self._nV)

    def getTransferFunction(self):
        return self._recTfunc

    def getVQTransferFunction(self, vCoeffs):
        #returns Vd*Qd
        vgTfunc=self._forTfunc*\
                self._dfilter.getTransferFunction(vCoeffs)
        return vgTfunc

    def getNvanCittert(self):
        return self._nV
#-----------------------------------------------------------------------
