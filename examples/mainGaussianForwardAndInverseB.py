from include import *

nVanCit=500
wavenumber=Wavenumber()
errors=Errors()

flt=Gaussian(wavenumber, param=4.0)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, adaptive=True, fixOne=[0.0], fixZero=[np.pi], error="msqe")
opt.run()

gdfor=dflt.getTransferFunction()
tfor=flt.getTransferFunction()
mpercF=errors.maperc(gdfor, tfor)
msqeF=errors.msqe(gdfor, tfor)
print("Forward filter:")
print("MSQE=", msqeF)
print("MPERC=", mpercF)

qFunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
vdflt=DiscreteFilter(wavenumber, isForward=False, nStencil=29)
iopt=Optimiser(qFunc, vdflt, adaptive=True, optError=1.0e-5, fixOne=[0.0], monotonePos=True, error="msqe")
iopt.run()

q=qFunc.getTransferFunction()
vg=vdflt.getTransferFunction()*gdfor
maxOverShootPerc=(np.max(vg)-1.0)*100.0
mpercI=errors.maperc(vg, q)
msqeI=errors.msqe(vg, q)
print("Inverse filter:")
print("MSQE=", msqeI)
print("MPERC=", mpercI)
print("Max over shoot % (VG)=", maxOverShootPerc)

out=Output(flt, dflt, qFunc, vdflt)
out.write()

plotAll(wavenumber.get(), flt, dflt, qFunc, vdflt, nVanCit)

