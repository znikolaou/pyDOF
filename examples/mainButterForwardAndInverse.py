from include import *

nVanCit=50
wavenumber=Wavenumber()
errors=Errors()

flt=Butterworth(wavenumber, paramKc=1.0197, paramN=6, fType="lowpass")
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, adaptive=True, monotoneNeg=True, fixOne=[0.0], fixZero=[np.pi], error="msqe")
opt.run()

gdfor=dflt.getTransferFunction()
tfor=flt.getTransferFunction()
mpercF=errors.maperc(gdfor, tfor)
msqeF=errors.msqe(gdfor, tfor)
print("Forward filter:")
print("MSQE=", msqeF)
print("MPERC=", mpercF)

qFunc=ReconstructedTransferFunction(dflt, nVanCittert=nVanCit)
mf=dflt.getNs()
vdflt=DiscreteFilter(wavenumber, isForward=False, nStencil=dflt.getNs())
iopt=Optimiser(qFunc, vdflt, adaptive=True, monotonePos=True, optError=1.0e-5, fixOne=[0.0], error="msqe")
iopt.run()

q=qFunc.getTransferFunction()
vg=vdflt.getTransferFunction()*gdfor
v=vdflt.getTransferFunction()
maxOverShootPerc=(np.max(vg)-1.0)*100.0
#mpercI=errors.maperc(vg, q)
#msqeI=errors.msqe(vg, q)
#mpercI=errors.maperc(v, q/g)
#msqeI=errors.msqe(vg, q/g)


#print("Inverse filter:")
#print("MSQE=", msqeI)
#print("MPERC=", mpercI)
#print("Max over shoot % (VG)=", maxOverShootPerc)

#out=Output(flt, dflt, qFunc, vdflt)
#out.write()

plotAll(wavenumber.get(), flt, dflt, qFunc, vdflt, nVanCit)

