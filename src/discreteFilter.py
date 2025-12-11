from utils import *
#-----------------------------------------------------------------------
class DiscreteFilter(Log, Const):
    """
    Class for representing a discrete filter on a given 
    stencil.

    Attributes: 
        wavenumber (Wavenumber): A wavenumber object.
        nStencil (int): The filter stencil size 2*nm+1 spanning [-nm, nm].
        isForward (bool): True for forward filter, False for inverse.

    Author: Z. Nikolaou.
    """

    def __init__(self, wavenumber,  nStencil=3, isForward=True):
        Log.__init__(self, type(self).__name__)
        self._wavenumber=wavenumber
        self._ns=nStencil
        self._isForward=isForward
        self.__checkSymmetric()
        self._k=self._wavenumber.get()
        self._nk=self._wavenumber.getNk()
        self.__checkSufficientStencil()
        self.coeffs=None
        self.optError=None

    def __checkSymmetric(self):
        if self._ns%2==0: 
            self.printError("__checkSymmetric", "stencil size not symmetric ")
            exit()

    def __checkSufficientStencil(self):
        if self._ns<3: 
            self.printError("_checkSufficientStencil", "stencil size<3")
            exit()

    def isForward(self): return self._isForward

    def getAllCoeffs(self):
        nm=len(self.coeffs)-1
        allCoeffs=np.zeros(2*nm+1)
        allCoeffs[nm:2*nm+1]=self.coeffs
        for i in range(0, nm): allCoeffs[i]=self.coeffs[nm-i]
        return allCoeffs

    def getK(self): return self._k

    def getNk(self): return self._nk

    def getNm(self): return int((self._ns-1)/2)

    def getNmPlusOne(self): return self.getNm()+1 

    def getNs(self): return self._ns

    def __setNsFromNmPlusOne(self, nmp1): self._ns=2*(nmp1-1)+1

    def setCoeffs(self, coeffs):
        """
        Sets coeffs g0, g1, ..., gm and total stencil 
        size ns accordingly. 
        """
        self.coeffs=coeffs
        self.__setNsFromNmPlusOne(len(self.coeffs))
    
    def setOptError(self, error): self.optError=error

    def getTransferFunction(self, halfCoeffs=None, nVanCittert=0):
        if(halfCoeffs is None): halfCoeffs=self.coeffs
        hg=np.zeros(self._nk)
        nh=len(halfCoeffs)
        for i in range(self._nk): 
            hg[i]=hg[i]+halfCoeffs[0]
            for j in range(1, nh):
                hg[i]=hg[i]+Const.TWO*halfCoeffs[j]*np.cos(j*self._k[i])

        return self.__getReconstructedTransferFunction(hg,nVanCittert)

    def __getReconstructedTransferFunction(self, tFunc, nVanCittert):
        if(nVanCittert==0): return tFunc
        q=np.zeros(self._nk)
        for i in range(self._nk):
            q[i]=self.ONE-(self.ONE-tFunc[i])**(nVanCittert+1)
        return q
#-----------------------------------------------------------------------
