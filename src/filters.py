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

    def getCutOff(self): return np.sqrt(24.0*np.log(Const.TWO))/self._param
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
class Inverse(Filter, Const):
    """
      Attributes: 

    Author: Z. Nikolaou
    """
 
    def __init__(self, forFilter, nV=0):
        self._forFilter=forFilter
        self._forTfunc=self._forFilter.getTransferFunction()
        self._k=self._forFilter._k
        self._nV=nV
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Inverse"

    def __buildTransferFunction(self):
        nk=len(self._forTfunc)
        tFunc=np.zeros(nk)
        for i in range(nk):
            tFunc[i]=(Const.ONE-(Const.ONE-self._forTfunc[i])**(self._nV+Const.ONE))/self._forTfunc[i]
        return tFunc

    def getNvanCittert(self): return self._nV
#-----------------------------------------------------------------------
class Sharp(Filter):
    """
    Class for designing sharp filters (lowpass, highpass, bandpass, bandstop). 

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        fType (str): "lowpass", "highpass", "bandpass", "bandstop". 
        kbounds (list): wavenumber boundaries. 

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumberOb, fType="lowpass", kbounds=[]):
        self._k=wavenumberOb.get()
        self._fType=fType
        self._kbounds=kbounds
        self._param=kbounds
        self._nkb=len(kbounds)
        self.__checkBounds()
        self.__checkSize()
        self.__checkIsMonotone()
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Sharp window"

    def __checkBounds(self):
        maxkb=np.max(self._kbounds)
        minkb=np.min(self._kbounds)
        maxk=np.max(self._k)
        mink=np.min(self._k)
        ier=0
        if minkb<mink or minkb>maxk: ier=1
        if maxkb<mink or maxkb>maxk: ier=2
        if ier!=0:
            print('*Sharp.__checkBounds() error:')
            print('\t kbounds outside wavenumber range! Terminating ...')
            exit()
        return
 
    def __checkSize(self):
        ier=0
        if self._fType=="lowpass" or self._fType=="highpass": 
            if self._nkb!=1: ier=1
        else:
            if (self._nkb%2)!=0: 
                ier=2
        if ier!=0:
            print('Sharp.__checkSize() error:')
            print('\t len(kbounds) must be 1 for fType=lowpass/highpass and len(kbounds)')
            print('\t len(kbounds) must be even for fType=bandpass/bandstop.')
            print('\t Terminating ...')
            exit()
        return

    def __checkIsMonotone(self):
        if self._nkb==1: return
        isMonotone=True
        for i in range(self._nkb-1):
            if self._kbounds[i+1]<=self._kbounds[i]: 
                isMonotone=False
                break
        if not isMonotone:
            print('*Sharp.__checkIsMonotone() error:')
            print('\t Supplied kbounds range should be monotonically increasing. Terminating ...')
            exit()
        return
                    
    def __buildTransferFunction(self):
        nk=len(self._k)
        if self._fType=="lowpass":
            tFunc=np.zeros(nk)
            for i in range(nk):
                if self._k[i]<=self._kbounds[0]: tFunc[i]=Const.ONE
        elif self._fType=="highpass":
            tFunc=np.ones(nk)
            for i in range(nk):
                if self._k[i]<=self._kbounds[0]: tFunc[i]=Const.ZERO
        else:
            indices=[self.__getFirst(kb, self._k) for kb in self._kbounds]
            ni=len(indices)
            if self._fType=="bandpass":
                tFunc=np.zeros(nk)
                for i in range(0, ni, 2): tFunc[indices[i]:indices[i+1]]=Const.ONE
            elif self._fType=="bandstop":
                tFunc=np.ones(nk)
                for i in range(0, ni, 2): tFunc[indices[i]:indices[i+1]]=Const.ZERO

        return list(tFunc)

    def __getFirst(self, x, xList):
        for i in range(len(xList)):
            if x<xList[i]:
                ip=i
                break
        return ip

    def getCutOff(self): return [self._param[0], self._param[1]]
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
    Class for Butterworth filter: ^G(k)=1/sqrt(1+3(k/kc)**2n)

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        paramN  (double) : sharpness 
        paramKc (double) : cut-off wavenumber

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumber, paramKc, paramN, fType=None):
        self._k=wavenumber.get()
        self._param=[paramKc, paramN]
        self._fType=fType
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Butterworth"

    def __buildTransferFunction(self):
        tFunc=[]
        exp=2.0e0*self._param[1]
        for k in self._k:
            hg=1.0e0/np.sqrt(1.0e0+3.0e0*(k/self._param[0])**exp)
            tFunc.append(hg)
        if self._fType=='highPass':
            return np.ones(len(self._k))-tFunc
        elif self._fType=='bandPass':
            pass
        else:
            return tFunc

    def getCutOff(self): return self._param[0]

    def getOrder(self): return self._param[1]

    def getKfromG(self, v): return self.getCutOff()*(1.0/3.0*(1.0/(v*v)-1.0))**(1./(2.0*self.getOrder()))

    def getKpass(self): return self.getKfromG(0.95)

    def getKstop(self): return self.getKfromG(0.05)

    def getTransWidth(self): return self.getKstop()-self.getKpass()

#-----------------------------------------------------------------------
class Constant(Filter):
    """
    Class for constant filter: ^G(k)=v

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object.
        v (float): filter constant value

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumber, v=1.0e0):
        self._k=wavenumber.get()
        self._param=v
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Constant"

    def __buildTransferFunction(self):
        tFunc=[]
        for k in self._k: tFunc.append(self._param)
        return tFunc
#-----------------------------------------------------------------------

class ErrorFunction(Filter):
    """
    Class for error-function-based filter: ^G(k)=0.5*(1-erf((kh-kch)sqrt(pi)/a))

    Attributes: 
        wavenumber (Wavenumber): the wavenumber object. 
        paramA  (double): sharpnesss
        paramKc (double): cut-off wavenumber

    Author: Z. Nikolaou
    """
    def __init__(self, wavenumber, paramKc, paramA):
        self._k=wavenumber.get()
        self._param=[paramKc, paramA]
        self._tFunc=self.__buildTransferFunction()
        Filter.__init__(self, self._k, self._tFunc)

    def name(self): return "Error-function"

    def __buildTransferFunction(self):
        tFunc=[]
        exp=2.0e0*self._param[1]
        for k in self._k:
            arg=(k-self._param[0])*np.sqrt(np.pi)/self._param[1]
            hg=0.5e0*(1.0e0-math.erf(arg))
            tFunc.append(hg)
        return tFunc

    def getCutOff(self): return self._param[0]
#-----------------------------------------------------------------------

