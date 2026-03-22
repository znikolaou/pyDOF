from include import *

nVanCit=50

wavenumber=Wavenumber()
errors=Errors()

flt=Butterworth(wavenumber, 1.01967, 4, 'highPass')
dflt=DiscreteFilter(wavenumber, nStencil=31)
opt=Optimiser(flt, dflt)
opt.run()
percErrorF=errors.maperc(dflt.getTransferFunction(), flt.getTransferFunction())
print("Forward filter mean abs perc error=", percErrorF)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5, isNoDampingAtK=[np.pi])
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
maxOverShootPerc=errors.maxOvershoot(vgTfunc)
percErrorI=errors.maperc(vgTfunc, recTfunc.getTransferFunction())
print("Inverse filter mean abs perc error=", percErrorI)
print("Max over shoot % (reconstructed)=", maxOverShootPerc)

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

