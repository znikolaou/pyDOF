from include import *

nVanCit=200

wavenumber=Wavenumber()
errors=Errors()

flt=Butterworth(wavenumber, np.pi/2.0, 10, 'lowPass')
dflt=DiscreteFilter(wavenumber, nStencil=55)
opt=Optimiser(flt, dflt, isNoDampingAtK=[0.0], isAdaptive=False, isPositive=True, isMonotonicNegative=True)
opt.run()
percErrorF=errors.maperc(dflt.getTransferFunction(), flt.getTransferFunction())
print("Forward filter mean abs perc error=", percErrorF)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, nStencil=151, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5, isNoDampingAtK=[0.0], isMonotonicPositive=True)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
maxOverShootPerc=errors.maxOvershoot(vgTfunc)
percErrorI=errors.maperc(vgTfunc, recTfunc.getTransferFunction())
print("Inverse filter mean abs perc error=", percErrorI)
print("Max over shoot % (reconstructed)=", maxOverShootPerc)

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

