from constants import *
from utils import *
import os
#-----------------------------------------------------------------------
class Output(Const):
    """
    Class for writing to file. 
    
    Attributes: 
        targetFilter (Filter) : Target filter object.
        discreteFilter (DiscreteFilter) : The discrete filter object.
        reconstructedFilter (ReconstructedTransferFunction) : Target reconstructed filter.
        discreteInverseFilter (DiscreteFilter) : The discrete filter object for the inverse filter.
        outputDir (string) : Output directory. 
        fileName (string)  : Name of output file. 

    Author: Z. Nikolaou. 
    """
    def __init__(self, targetFilter, discreteFilter, reconstructedFilter=None, discreteInverseFilter=None, \
            outputDir="./output/", fileName="coeffs"):
        self._dir=outputDir
        self._fileName=fileName
        self.__createDir()
        self._target=targetFilter
        self._dfilter=discreteFilter
        self._rfilter=reconstructedFilter
        self._ifilter=discreteInverseFilter
        self._fullPath=self._dir+"/"+self._fileName+".txt"

    def __createDir(self):
        isExists=os.path.exists(self._dir)
        if not isExists: os.makedirs(self._dir)
        return

    def write(self):
        f=open(self._fullPath,'w')
        coeffs=self._dfilter.coeffs
        nmp1=len(coeffs)
        nm=nmp1-1
        f.write(self.HEADER_OUT)
        f.write('\n#Forward filter')
        f.write('\n! Name='+self._target.name())
        nprm=np.size(self._target._param)
        if nprm==1: 
            f.write('\n! Parameter='+Form.short(self._target._param))
        else:
            for i in range(nprm): 
                f.write('\n! Parameter='+Form.short(self._target._param[i]))
        #f.write('\n! k_ch=' \
        #        +Form.short(self._target.getCutOff()))
        #f.write('\n! Damping at pi='+ \
        #        Form.short(self._target.getNyquistDamping()))
        f.write('\n! Stencil size='+str(self._dfilter.getNs()))
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

        return

    def __writeInverse(self, f):
        coeffs=self._ifilter.coeffs
        nmp1=len(coeffs)
        f.write('\n#Inverse filter')
        f.write('\n! Reconstruction order (N van Cittert)='+ \
                str(self._rfilter.getNvanCittert()))
        f.write('\n! Stencil size='+str(self._ifilter.getNs()))
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
