from include import *

nVanCit=5

wavenumber=Wavenumber()
errors=Errors()

flt=Gaussian(wavenumber, param=4.0)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, isAdaptive=True, isNoDampingAtZero=True, isMonotonicNegative=True)
opt.run()
percErrorF=errors.maperc(dflt.getTransferFunction(), flt.getTransferFunction())
print("Forward filter mean abs perc error=", percErrorF)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5, isNoDampingAtZero=True, isMonotonicPositive=True)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
maxOverShootPerc=errors.maxOvershoot(vgTfunc)
percErrorI=errors.maperc(vgTfunc, recTfunc.getTransferFunction())
print("Inverse filter mean abs perc error=", percErrorI)
print("Max over shoot % (reconstructed)=", maxOverShootPerc)

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

