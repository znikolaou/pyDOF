from utils import *
import datetime
class FortranGenerator(Const):
    """
    Class for generating a Fortran source code containing info 
    on the discrete filter(s).  

    Attributes: 
        path (string): 
            Path to store header file.
        discreteFiltersList (Discrete Filter []):
            List of DiscreteFilter objects.

    Author: Z. Nikolaou.
    """

    # TODO: CLEAN-UP 

    def __init__(self, path, discreteFiltersList):
        self._path=path
        self._filtersList=discreteFiltersList
        self._nFilters=len(self._filtersList)
        self._ncList=self.__getNcList()
        self._ncm1List=[i-1 for i in self._ncList]
        self._isForwardList=self.__getIsForwardList()
        self.__fl=open(path+'/'+self.FORTRAN_HEADER,'w')
    
    def __getNcList(self):
        return [len(ob.coeffs) for ob in self._filtersList]

    def __getIsForwardList(self):
        return [1 if ob.isForward() else 0 for ob in self._filtersList]

    def generateHeader(self):
        NC='NC'
        NCM1='NCM1'
        NCMAX='NC_MAX'
        NCMAXM1='NC_MAX_M1'
        NFILTERS='NFILTERS'
        COEFS='COEFS'

        ncMax=max(self._ncList)
        ncMaxM1=ncMax-1

        self.__writeComment('---------------------------')
        self.__writeComment('This code has been auto-generated:')
        self.__writeComment(str(datetime.datetime.now()))
        self.__writeComment('')
        self.__writeComment(self.LOGO)
        self.__writeComment('Author: Z. Nikolaou')
        self.__writeComment('---------------------------')
        self.__writeLine('INTEGER '+NFILTERS+','+NCMAX+','+NCMAXM1)
        self.__writeLine('PARAMETER('+NFILTERS+'='+str(self._nFilters) \
                                 +','+NCMAX+'='+str(ncMax) \
                                 +','+NCMAXM1+'='+str(ncMaxM1)+')')
        self.__writeLine('INTEGER '+NC+'('+NFILTERS+')')
        self.__writeLine('INTEGER '+NCM1+'('+NFILTERS+')')
        self.__writeLine('INTEGER ISFORWARD('+NFILTERS+')')
        self.__writeLine('DOUBLE PRECISION '+COEFS+'(0:'+NCMAXM1+','+ \
                NFILTERS+')')
        self.__writeData(NC,self._ncList)
        self.__writeData(NCM1,self._ncm1List)
        self.__writeData('ISFORWARD',self._isForwardList)
        #Filter coefs
        nFilters=len(self._filtersList)
        coeffsArray=np.zeros(ncMax)
        for j in range(nFilters):
            ob=self._filtersList[j]
            nc=self._ncList[j]
            coeffsArray[:]=0.0e0
            coeffsArray[0:nc]=ob.coeffs[0:nc]
            self.__writeData(COEFS+'(:,'+str(j+1)+')', coeffsArray, \
                    self.NUM_FORMAT)
        self.__writeComment('------END------')
        self.__fl.close()

    def __writeData(self, string, varList, fmt=None):
        self.__writeLine('DATA '+string+'/')
        nm1=len(varList)-1
        for i in range(nm1):
            self.__writeContinuationLine(self.__getFormattedString(varList[i], fmt)+',')
        self.__writeContinuationLine(self.__getFormattedString(varList[nm1], fmt)+'/')

    def __getFormattedString(self, value, fmt=None):
        return str(value) if fmt is None else str(fmt(value))
        
    def __getLine(self, line):
        return '      '+line+'\n'

    def __getCommentLine(self, line):
        return '!     '+line+'\n'

    def __getContinuationLine(self, line):
        return '     & '+line+'\n'

    def __writeLine(self, line):
        self.__fl.write(self.__getLine(line))

    def __writeComment(self, line):
        self.__fl.write(self.__getCommentLine(line))

    def __writeContinuationLine(self, line):
        self.__fl.write(self.__getContinuationLine(line))






