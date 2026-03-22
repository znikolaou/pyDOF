from include import *

nVanCit=500

wavenumber=Wavenumber()
errors=Errors()

flt=Butterworth(wavenumber, np.pi/4.0, 4, 'lowPass')
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, isNoDampingAtK=[0.0], isPositive=True, isMonotonicNegative=True)
opt.run()
percErrorF=errors.maperc(dflt.getTransferFunction(), flt.getTransferFunction())
print("Forward filter mean abs perc error=", percErrorF)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, nStencil=dflt.getNs(), isForward=False)
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

