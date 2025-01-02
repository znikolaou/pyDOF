from constants import *
from utils import *
#-----------------------------------------------------------------------
class Output(Const):
    """
    Class for writing to file. 
    
    Attributes: 
        outputDir: string 
            Output directory. 
        fileName: string
            Name of output file. 
        coeffs: double [] 
            Discrete filter coeffs g0, g1, ..., gm.
        targetFilter: Filter 
            Target filter.
        discreteFilter: DiscreteFilter
            The discrete filter object.
        reconstructedFilter: ReconstructedTransferFunction 
            Target reconstructed filter.
        inverseFilter: DiscreteFilter
            The discrete filter object for the inverse filter.

    Author: Z. Nikolaou. 
    """
    def __init__(self, outputDir, fileName, targetFilter, \
            discreteFilter, reconstructedFilter=None, \
            discreteInverseFilter=None):
        self._dir=outputDir
        self._fileName=fileName
        self._target=targetFilter
        self._dfilter=discreteFilter
        self._rfilter=reconstructedFilter
        self._ifilter=discreteInverseFilter
        self._fullPath=self._dir+"/"+self._fileName+".txt"

    def write(self):
        f=open(self._fullPath,'w')
        coeffs=self._dfilter.coeffs
        nmp1=len(coeffs)
        nm=nmp1-1
        f.write(self.HEADER_OUT)
        f.write('\n#Forward filter')
        f.write('\n! Name='+self._target.name())
        f.write('\n! Parameter='+Form.short(self._target._param))
        f.write('\n! k_fh (\hat{G}(k_fh)=0.5)=' \
                +Form.short(self._target.getCutOff()))
        f.write('\n! Damping at pi (\hat{G}(\pi)='+ \
                Form.short(self._target.getNyquistDamping()))
        f.write('\n! Stencil size (2*nm_f+1)='+str(self._dfilter.getNs()))
        f.write('\n! Optimisation error='+ \
                Form.short(self._dfilter.optError))
        f.write('\n! Coeffs g[0], g[1], ... , g[nm_f]:')
        f.write('\n'+str(nmp1))
        for j in range(nmp1): \
                f.write('\n'+str(j)+" "+Form.long(coeffs[j]))
        f.write('\n'+self.LINE)
        if(self._rfilter != None):
            if(self._ifilter != None): 
                self.__writeInverse(f)
        f.close()

    def __writeInverse(self, f):
        coeffs=self._ifilter.coeffs
        nmp1=len(coeffs)
        f.write('\n#Inverse filter')
        f.write('\n! van Cittert iterations='+ \
                str(self._rfilter.getNvanCittert()))
        f.write('\n! Stencil size (2*nm_i+1)='+str(self._ifilter.getNs()))
        f.write('\n! Optimisation error='+ \
                Form.short(self._ifilter.optError))
        f.write('\n! Coeffs g[0], g[1], ... , g[nm_i]:')
        f.write('\n'+str(nmp1))
        for j in range(nmp1): \
                f.write('\n'+str(j)+" "+Form.long(coeffs[j]))
        f.write('\n'+self.LINE)
        f.close()

        return
#-----------------------------------------------------------------------
