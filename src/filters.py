from constants import *
from utils import *
from scipy.optimize import fsolve
#-----------------------------------------------------------------------
class Filter(Log, Const):
    """
    Base class for representing continuous filters. 

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
      
    def getNyquistDamping(self): return self._tFunc[len(self._k)-1]

    def getTotalDamping(self): return np.trapz(self._tFunc, self._k)

    def getDampingEfficiency(self):
        lk=max(self._k)-min(self._k)
        return self.getDamping()/lk

    @abstractmethod 
    def name(self): pass

    @abstractmethod
    def getNominalStencil(self): pass

    @abstractmethod
    def get(self): pass

#-----------------------------------------------------------------------
class Gaussian(Filter):
    """
    Class for Gaussian filter: ^G(k)=exp(-Delta**2k**2/24)

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        param (double): Delta/h ratio for controlling the filter shape. 

    Author: Z. Nikolaou
    """
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
    """
    Class for the Helmholtz filter: ^G(k)=2/(kL)sin(kL/2)

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        param (double): lambda for controlling the filter shape. 

    Author: Z. Nikolaou
    """
 
    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Helmholtz"

    def getNominalStenicl(self): return int(2.0*7.0*self._param+1)

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            hg=Const.ONE/(Const.ONE+self._param**Const.TWO*k**Const.TWO)
            tFunc.append(hg)
        return tFunc

    def getCutOff(self): return Const.ONE/self._param
#-----------------------------------------------------------------------
class Implicit(Filter):
    """
    Class for a 3-point implicit filter: ^G(k)=(0.5+a)(1+cos(kh))/(1+2acos(kh))

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        param (double): parameter for controlling the filter shape. 

    Author: Z. Nikolaou
    """
 
    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "3-point Implicit"

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            num=(Const.HALF+self._param)*(Const.ONE+np.cos(k))
            den=Const.ONE+Const.TWO*self._param*np.cos(k)
            tFunc.append(num/den)
        return tFunc

    def getCutOff(self): return np.arccos(-Const.TWO*self._param)
#-----------------------------------------------------------------------
class SharpCutOff(Filter):
    """
    Class for a shartp cut-off filter.

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        param (double): filter cut-off wavenumber (kch). 

    Author: Z. Nikolaou
    """
 
    def __init__(self, wavenumberOb, param):
        self._k=wavenumberOb.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Sharp cut-off"

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            if(k<=self._param):
                tFunc.append(Const.ONE)
            else:
                tFunc.append(Const.ZERO)
        return tFunc

    def getCutOff(self): return self._param

    def get(self, x):
        n=len(x)
        g=np.zeros(n)
        for i in range(n):
            if x[i] != 0.0e0:
                g[i]=1.0e0/np.pi*np.sin(self._param*x[i])/x[i]
            else:
                x[i]=self._param/np.pi
        return g
#-----------------------------------------------------------------------
class Box(Filter):
    """
    Class for Box filter: ^G(kh)=2/(kh L/h)sin(kh (L/h)/2)

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        param (double): L/h the region over which G(x)=1/L. 

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumber, param):
        self._k=wavenumber.get()
        self._param=param
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Box"

    def getNominalStencil(self): return int(self._param)+1

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k:
            if k!=0.0e0:
                hg=2.0e0/(k*self._param)*np.sin(k*self._param/2.0e0)
            else:
                hg=1.0e0
            tFunc.append(hg)
        return tFunc

    def getCutOff(self): return fsolve(self.funcForCutOff, 0.0)

    def funcForCutOff(self, x): return np.sin(x*self._param/2.0e0)-x/4.0e0*self._param

    @staticmethod 
    def getParamForBoxFilterCutOff(kc): return fsolve(Box.funcForBoxFilterParam, np.pi, args=kc)

    @staticmethod
    def funcForBoxFilterParam(x, kc): return np.sin(kc*x/2.0e0)-kc/4.0e0*x

#-----------------------------------------------------------------------
class Butterworth(Filter):
    """
    Class for Butterworth filter: ^G(k)=1/sqrt(1+(k/kc)**2n)

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        paramN  (double)
        paramKc (double)

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumber, paramKc, paramN):
        self._k=wavenumber.get()
        self._param=[paramKc, paramN]
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Butterworth"

    def __buildTransferFunction(self):
        tFunc=[]
        exp=2.0e0*self._param[1]
        for k in self._k:
            hg=1.0e0/np.sqrt(1.0e0+3.0e0*(k/self._param[0])**exp)
            tFunc.append(hg)
        return tFunc

    def getCutOff(self): return self._param[0]
#-----------------------------------------------------------------------

