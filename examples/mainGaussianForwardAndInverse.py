from include import *

nVanCit=5

wavenumber=Wavenumber()
errors=Errors()

flt=Gaussian(wavenumber, param=4.0)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, isAdaptive=True, isNoDampingAtZero=True, isFullDampingAtPi=True, isMonotonicNegative=True)
opt.run()
print("Forward filter perc error=", percError)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5, isNoDampingAtZero=True, isMonotonicPositive=True)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
maxOverShootPerc=errors.maxOvershoot(vgTfunc)
percErrorF=errors.mperc(flt.getTransferFuction(), dflt.getTransferFunction())
print("Max over shoot % (reconstructed)=", maxOverShoot, maxOverShootPerc)

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

