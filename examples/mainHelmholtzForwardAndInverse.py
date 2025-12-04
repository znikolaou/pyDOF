from include import *

wavenumber=Wavenumber()

flt=Helmholtz(wavenumber, param=1.0/1.01967)
dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, isZeroDampingAtPi=False)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=300)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, isZeroDampingAtPi=False, optError=1.0e-5)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}$")
plt.plot(k, dflt.getTransferFunction(), 'r', label="$\hat{G}_d$")
plt.plot(k, recTfunc.getTransferFunction(), 'bd', mfc="None", label="${\hat{Q}_d}^N$")
plt.plot(k, vgTfunc, 'm', label="$\hat{V}_d\hat{G}_d$")
plt.plot(k, vdflt.getTransferFunction()/np.max(vdflt.getTransferFunction()), 'g', label="$\hat{V}_d+$")
plt.xlabel('$kh$')
plt.legend(frameon=False)
plt.grid('on')
plt.show()


