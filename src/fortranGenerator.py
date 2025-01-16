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
        self._PRM='PARAMETER'

    def __getNcList(self):
        return [len(ob.coeffs) for ob in self._filtersList]

    def __getIsForwardList(self):
        return [1 if ob.isForward() else 0 for ob in self._filtersList]

    def generateHeader(self):
        NCOEF='NCOEF'
        NCOEFM1='NCOEFM1'
        NCMX='NCMX'
        NGMX='NGMX'
        NF='NF'
        COEFS='COEFS'
        NXMX='NXMX'
        NYMX='NYMX'
        NZMX='NZMX'
        NCM1OR1='NCM1OR1'
        nxmx=1000

        ncMax=max(self._ncList)
        ncMaxM1=ncMax-1

        self.__writeComment('----------------------------------')
        self.__writeComment('This code has been auto-generated:')
        self.__writeComment(self.__buildDate())
        self.__writeComment('')
        self.__writeComment(self.LOGO)
        self.__writeComment('Author: Z. Nikolaou')
        self.__writeComment('----------------------------------')

        self.__writeIntegers([NF, NCMX, NGMX])
        self.__writeParam(NF, str(self._nFilters))
        self.__writeParam(NCMX, str(ncMax))
        self.__writeParam(NGMX, str(ncMaxM1))
        self.__writeIntegerArray(NCOEF, NF)
        self.__writeIntegerArray(NCOEFM1, NF)
        self.__writeIntegerArray('ISFORWARD', NF)
        #todo:create double array method
        self.__writeLine('DOUBLE PRECISION '+COEFS+'(0:'+NGMX+','+ \
                NF+')')

        self.__writeData(NCOEF,self._ncList)
        self.__writeData(NCOEFM1,self._ncm1List)
        self.__writeData('ISFORWARD',self._isForwardList)

        #todo: reformat/clean
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
        self.__writeIntegers([NXMX,NYMX,NZMX])
        self.__writeParam(NXMX, str(nxmx))
        self.__writeParam(NYMX, str(nxmx))
        self.__writeParam(NZMX, str(nxmx))
        self.__writeLine('INTEGER '+NCM1OR1+'('+str(3)+','+NF+')')
        self.__writeLine('INTEGER '+'IL(-'+NGMX+':'+NGMX+','+NXMX+','+NF+')')
        self.__writeLine('INTEGER '+'JL(-'+NGMX+':'+NGMX+','+NYMX+','+NF+')')
        self.__writeLine('INTEGER '+'KL(-'+NGMX+':'+NGMX+','+NZMX+','+NF+')')       
        self.__writeLine('DOUBLE PRECISION '+'COEFS_A(-'+NGMX+':'+NGMX+','+ \
                                                     '-'+NGMX+':'+NGMX+','+ \
                                                     '-'+NGMX+':'+NGMX+',NF'+')')
        self.__writeLine('COMMON/PYDOF_COMM/ IL,JL,KL,COEFS_A,'+NCM1OR1)
        self.__writeComment('------END------')
        #close file
        self.__fl.close()

    def __buildDate(self):
        return str(datetime.datetime.now())

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

    def __writeParam(self, name, size):
        return self.__writeLine(self._PRM+' ('+name+'='+size+')')

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






