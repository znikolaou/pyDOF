from utils import *
import datetime
class FortranGenerator(Const):
    """
    Class for generating a Fortran header file
    containing the discrete filter coeffs. 

    Attributes: 
        path (string): path to store header file
        filtersList (Discrete Filter []): list of DiscreteFilter objects
        forwardOrInverseList (int []): integer list specifying 
                                       forward (1) or inverse(0) filter

    Author: Z. Nikolaou.
    """
    def __init__(self, path, filtersList, forwardOrInverseList):
        self._path=path
        self._filtersList=filtersList
        self._fORi=forwardOrInverseList
        self._nFilters=len(filtersList)
        self._ncMax=self.__getNcMax()
        self.__fl=open(path+'/'+'pydof_common.h','w')
    
    def __getNcMax(self):
        return max([len(ob.coeffs) for ob in self._filtersList])

    def generateHeader(self):
        NC='NC'
        NCMAX='NC_MAX'
        NFILTERS='NFILTERS'
        COEFS='COEFS'
        NCMAXM1=str(self._ncMax-1)
        self.__writeComment('---------------------------')
        self.__writeComment('This code has been auto-generated:')
        self.__writeComment(str(datetime.datetime.now()))
        self.__writeComment('')
        self.__writeComment(self.LOGO)
        self.__writeComment('Author: Z. Nikolaou')
        self.__writeComment('---------------------------')
        self.__writeLine('INTEGER '+NFILTERS+','+NCMAX)
        self.__writeLine('PARAMETER('+NFILTERS+'='+str(self._nFilters) \
                +','+NCMAX+'='+str(self._ncMax)+')')
        self.__writeLine('INTEGER '+NC+'('+NFILTERS+')')
        self.__writeLine('INTEGER ISFORWARD('+NFILTERS+')')
        self.__writeLine('DOUBLE PRECISION '+COEFS+'(0:'+NCMAXM1+','+NFILTERS+')')
        

        nmAll=[len(ob.coeffs) for ob in self._filtersList]
        self.__writeParameters(NC,nmAll)

        self.__writeParameters('ISFORWARD',self._fORi)
        
        nFilters=len(self._filtersList)
        for j in range(nFilters):
            self.__writeLine('PARAMETER(')
            ob=self._filtersList[j]
            nm=len(ob.coeffs)
            for i in range(nm-1):
                coeff=Const.NUM_FORMAT(ob.coeffs[i])
                self.__writeContinuationLine(COEFS+\
                        '('+str(i)+','+str(j+1)+')='+str(coeff)+',')
            i=i+1
            coeff=Const.NUM_FORMAT(ob.coeffs[nm-1])
            self.__writeContinuationLine(COEFS+'('+str(i)+',' \
                    +str(j+1)+')='+str(coeff)+')')

        self.__writeComment('------END------')
        self.__fl.close()

    def __writeParameters(self, string, varList):
        self.__writeLine('PARAMETER(')
        k=1
        n=len(varList)
        nm1=n-1
        for i in range(nm1):
            self.__writeContinuationLine(string+'('+str(k)+')='+str(varList[i])+',')
            k=k+1
        self.__writeContinuationLine(string+'('+str(k)+')='+str(varList[nm1])+')')

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






