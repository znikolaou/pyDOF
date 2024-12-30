class FortranCodeGenerator(Const):

    _headerName="pydof_filterCommon.h"
    _sixSpace="      "

    def __init__(self, path, dFilterOb):
        self._path=path
        self._dFilterOb=dFilterOb
        self._fileName=self._path+"/"+self._headerName
        self._file=open(self._file, "w")
        self._generateCode()
        self._file.close()

    def _generateCode(self):
        self.writeLine("!THIS IS AUTO-GENERATED CODE")
        self.writeLine("!"+self.LOGO)
        self._writeLine("INTEGER NXM="+str(self._dfilterOb.getNm()))
        self._writeLine("DOUBLE PRECISION COEFFS(-NXM:NXM)")
        for i in range(self._dFilterOb.getNm()):
            self._writeLine("COEFFS("+str(i)+")="+ \
                    str(self._dFilteOb.coeffs[i]))
        self._writeLine("COMMON/COEFFS_BLK/ NXM,COEFFS")

    def _writeLine(self, line):
        self._file.write(self._sixSpace+line)






