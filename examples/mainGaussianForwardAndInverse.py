from include import *

wavenumber=Wavenumber()

nVanCit=500

flt=Gaussian(wavenumber, param=4.0)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-5)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
maxOverShoot=np.max(vgTfunc)
maxOverShootPerc=(maxOverShoot-1.0e0)*100.0
print("Max over shoot (reconstructed)=", maxOverShoot, "(", maxOverShootPerc,"%)")

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, recTfunc, vdflt, nVanCit)

