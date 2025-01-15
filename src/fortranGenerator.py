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
        self._INT='INTEGER'
        self._DBL='DOUBLE PRECISION'
        self._DATA='DATA'
    
    def __getNcList(self):
        return [len(ob.coeffs) for ob in self._filtersList]

    def __getIsForwardList(self):
        return [1 if ob.isForward() else 0 for ob in self._filtersList]

    def generateHeader(self):
        NC='NCOEF'
        NCM1='NCOEFM1'
        NCMAX='NCMX'
        NCMAXM1='NGMX'
        NFILTERS='NF'
        COEFS='COEFS'
        NXMAX='1000000'
        NYMAX=NXMAX
        NZMAX=NXMAX

        ncMax=max(self._ncList)
        ncMaxM1=ncMax-1

        self.__writeComment('----------------------------------')
        self.__writeComment('This code has been auto-generated:')
        self.__writeComment(str(datetime.datetime.now()))
        self.__writeComment('')
        self.__writeComment(self.LOGO)
        self.__writeComment('Author: Z. Nikolaou')
        self.__writeComment('----------------------------------')

        self.__writeIntegers([NFILTERS, NCMAX, NCMAXM1])
        self.__writeLine('PARAMETER('+NFILTERS+'='+str(self._nFilters) \
                                 +','+NCMAX+'='+str(ncMax) \
                                 +','+NCMAXM1+'='+str(ncMaxM1)+')')
        self.__writeIntegerArray(NC, NFILTERS)
        self.__writeIntegerArray(NCM1, NFILTERS)
        self.__writeIntegerArray('ISFORWARD', NFILTERS)
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
        #common block work arrays
        self.__writeComment('')
        self.__writeComment('COMMON BLOCK DATA')
        self.__writeComment('')
        self.__writeLine('INTEGER '+'IL(-'+NCMAXM1+':'+NCMAXM1+','+NXMAX+')')
        self.__writeLine('INTEGER '+'JL(-'+NCMAXM1+':'+NCMAXM1+','+NZMAX+')')
        self.__writeLine('INTEGER '+'KL(-'+NCMAXM1+':'+NCMAXM1+','+NYMAX+')')       
        self.__writeLine('DOUBLE PRECISION '+'COEFS_A(-'+NCMAXM1+':'+NCMAXM1+','+ \
                                                     '-'+NCMAXM1+':'+NCMAXM1+','+ \
                                                     '-'+NCMAXM1+':'+NCMAXM1+')')
        self.__writeLine('COMMON/PYDOF_COMM/ IL,JL,KL,COEFS_A')
        self.__writeComment('------END------')
        #close file
        self.__fl.close()

    def __writeData(self, string, varList, fmt=None):
        self.__writeLine(self._DATA+' '+string+'/')
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

    def __writeIntegers(self, namesList):
        self.__writeLine(self._INT+ \
                self.__buildCommaSeparated(namesList))

    def __writeDoubles(self, namesList):
        self.__writeLine(self._DBL+ \
                self.__buildCommaSeparated(namesList))
    
    def __writeIntegerArray(self, name, size):
        self.__writeLine(self._INT+self.__wrapInPars(name, size))

    def __writeDoubleArray(self, name, size):
        self.__writeLine(self._DBL+self.__wrapInPars(name, size))

    def __buildCommaSeparated(self, stringList):
        string=' '+stringList[0]
        n=len(stringList)
        for i in range(1, n):
            string=string+','+stringList[i]        
        return string

    def __wrapInPars(self, name, size):
        return ' '+name+'('+size+')'
     
    def __writeLine(self, line):
        self.__fl.write(self.__getLine(line))

    def __writeComment(self, line):
        self.__fl.write(self.__getCommentLine(line))

    def __writeContinuationLine(self, line):
        self.__fl.write(self.__getContinuationLine(line))






