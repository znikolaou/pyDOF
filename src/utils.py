from constants import *
from libs import *
#-----------------------------------------------------------------------
class Log():
    """
    Logger class.

    Attributes: 
        name (string): The name of the class/method for which to output
                       information.

    Author: Z. Nikolaou.
    """

    def __init__(self, name):
        self.__err="Error: "
        self.__warn="Warning: "
        self.__dot="."
        self.__nLine="\n"
        self.__name=name

    def __buildPre(self, mode):
        return self.__nLine+mode+self.__name+self.__dot

    def printError(self, method, message):
        print(self.__buildPre(self.__err)+method+" "+message)

    def printWarning(self, method, message):
        print(self.__buildPre(self.__warn)+method+" "+message)
#-----------------------------------------------------------------------
class Wavenumber(Log, Const):
    """
    Class for defining the wavenumber range. 

    Attributes:
        nk (int): the number of discretisation points for the 
                  wavenumber k spanning [0, pi].

    Author: Z. Nikolaou.
    """

    def __init__(self, nk=Const.NK_NOMINAL):
        Log.__init__(self, type(self).__name__) 
        self._nk=nk
        self._k=[]
        self.__set()

    def __set(self): self._k=np.linspace(self.ZERO, self.PI, self._nk)

    def get(self): return self._k

    def getNk(self): return self._nk
#-----------------------------------------------------------------------
class Form(Const):

    @staticmethod
    def short(number): return Const.NUM_FORMAT_SHORT(number)

    @staticmethod
    def long(number): return Const.NUM_FORMAT_LONG(number)
#-----------------------------------------------------------------------

class FUtil():

    @staticmethod
    def getCutOffIndex(g): return g.tolist().index(0.5)

    @staticmethod
    def getButterOrderForDeltaK(kc, deltaK, gpass=0.99, gstop=0.01): 
        sol=scipy.optimize.root(FUtil.funcButter, 10, args=(deltaK, kc, gpass, gstop))
        return sol.x

    @staticmethod
    def funcButter(n, deltaK, kc, gpass, gstop):
        kpass=kc*FUtil.getButterVal(n, gpass)
        kstop=kc*FUtil.getButterVal(n, gstop)
        return (kstop-kpass)-deltaK

    @staticmethod
    def getButterVal(n, g):
        return 1.0/(2.0*n)*np.log(1/3.0*(1/g**2.0-1.0))
    
        





