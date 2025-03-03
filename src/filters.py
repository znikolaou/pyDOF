from constants import *
from utils import *
#-----------------------------------------------------------------------
class Filter(Log, Const):
    """
    Parent class for representing continuous filters. 

    *Any new filter should inherit Filter! 

    Attributes: 
        k (double [])    : grid wavenumbers typically in [0, pi].
        tFunc (double []): corresponding transfer function.

    Author: Z. Nikolaou.
    """
    def __init__(self, k, tFunc):
        Log.__init__(self, type(self).__name__)
        self._k=k
        self._tFunc=tFunc
        self.__checkNyquistDamping()
    
    def __checkNyquistDamping(self):
        damping=self.getNyquistDamping() 
        if(abs(damping) > Const.TOL_NYQDAMP):
            self.printWarning("_checkNyqDamping",  \
                    "insufficient Nyquist damping "+str(damping))

    def getTransferFunction(self, nv=0):
        nk=len(self._k)
        tf=np.zeros(nk)
        for i in range(nk):
            tf[i]=Const.ONE-(Const.ONE-self._tFunc[i])**(nv+1)
        return tf
      
    def getNyquistDamping(self): 
        return self._tFunc[len(self._k)-1]

    def getTotalDamping(self):
        return np.trapz(self._tFunc, self._k)

    def getDampingEfficiency(self):
        lk=max(self._k)-min(self._k)
        return self.getDamping()/lk

    @abstractmethod 
    def name(self): pass

    @abstractmethod
    def getNominalStencil(self):
        pass
#-----------------------------------------------------------------------
class Gaussian(Filter):
    
    def __init__(self, wavenumber, param):
        self._k=wavenumber.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Gaussian"

    def getNominalStencil(self): return int(2.0*self._param+1)

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            hg=np.exp(-self._param**2.0*k**Const.TWO/24.0)
            tFunc.append(hg)
        return tFunc

    def getCutOff(self):
        return math.sqrt(24.0*math.log(Const.TWO))/self._param
#-----------------------------------------------------------------------
class Helmholtz(Filter):
 
    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self):
        return "Helmholtz"

    def getNominalStenicl(self): return int(2.0*7.0*self._param+1)

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            hg=Const.ONE/(Const.ONE+self._param**Const.TWO*k**Const.TWO)
            tFunc.append(hg)
        return tFunc

    def getCutOff(self):
        return Const.ONE/self._param
#-----------------------------------------------------------------------
class Implicit(Filter):

    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self):
        return "3-point Implicit"

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            num=(Const.HALF+self._param)*(Const.ONE+np.cos(k))
            den=Const.ONE+Const.TWO*self._param*np.cos(k)
            tFunc.append(num/den)
        return tFunc

    def getCutOff(self):
        return np.arccos(-Const.TWO*self._param)
#-----------------------------------------------------------------------
class SharpCutOff(Filter):

    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self):
        return "Sharp cut-off"

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            if(k<=self._param):
                tFunc.append(Const.ONE)
            else:
                tFunc.append(Const.ZERO)
        return tFunc
#-----------------------------------------------------------------------

