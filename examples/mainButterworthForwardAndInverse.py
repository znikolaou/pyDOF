from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

wavenumber=Wavenumber()
flt=Butterworth(wavenumber, 1.01967, 6.0)

dflt=DiscreteFilter(wavenumber, nStencil=55)
opt=Optimiser(flt, dflt, isFullDampingAtPi=False, optError=1.0e-6)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=300)
vdflt=DiscreteFilter(wavenumber, isForward=False, nStencil=89)
iopt=Optimiser(recTfunc, vdflt, isFullDampingAtPi=False, isNoDampingAtZero=False, isAdaptive=False, lamda=5.0)
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

