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
