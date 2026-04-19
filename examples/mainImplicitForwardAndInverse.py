from include import *

nVanCit=500

wavenumber=Wavenumber()
errors=Errors()

flt=Implicit(wavenumber, param=-0.262)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, adaptive=True, fixOne=[0.0], fixZero=[np.pi], monotoneNeg=True)
opt.run()
percErrorF=errors.maperc(dflt.getTransferFunction(), flt.getTransferFunction())
print("Forward filter mean abs perc error=", percErrorF)

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5, fixOne=[0.0], monotonePos=True)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
#maxOverShootPerc=errors.maxOvershoot(vgTfunc)
percErrorI=errors.maperc(vgTfunc, recTfunc.getTransferFunction())
print("Inverse filter mean abs perc error=", percErrorI)
#print("Max over shoot % (reconstructed)=", maxOverShootPerc)

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

