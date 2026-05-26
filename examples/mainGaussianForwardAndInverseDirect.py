from include import *

nVanCit=50
wavenumber=Wavenumber()
errors=Errors()

flt=Gaussian(wavenumber, param=6.0)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, adaptive=True, fixOne=[0.0], error="msqe")
opt.run()

gdfor=dflt.getTransferFunction()
tfor=flt.getTransferFunction()
mpercF=errors.maperc(gdfor, tfor)
msqeF=errors.msqe(gdfor, tfor)
print("Forward filter:")
print("MSQE=", msqeF)
print("MPERC=", mpercF)

targetInv=Inverse(dflt, nVanCit)
vdflt=DiscreteFilter(wavenumber)
iopt=Optimiser(targetInv, vdflt, adaptive=True, optError=1.0e-3, fixOne=[0.0], monotonePos=True, error="msqe")
iopt.run()

t=targetInv.getTransferFunction()
v=vdflt.getTransferFunction()
maxOverShootPerc=(np.max(v*dflt.getTransferFunction())-1.0)*100.0
#mpercI=errors.maperc(v, q)
#msqeI=errors.msqe(vg, q)
print("Inverse filter:")
#print("MSQE=", msqeI)
#print("MPERC=", mpercI)
print("Max over shoot % (VG)=", maxOverShootPerc)

#out=Output(flt, dflt, targetInv, vdflt)
#out.write()

plotAll(wavenumber.get(), flt, dflt, targetInv, vdflt, nVanCit)

