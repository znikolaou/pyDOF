from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

wavenumber=Wavenumber()
flt=Gaussian(wavenumber, 4.0)
nStencilNominal=flt.getNominalStencil()

dflt=DiscreteFilter(wavenumber, nStencil=nStencilNominal)
opt=Optimiser(flt, dflt)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=100)
vdflt=DiscreteFilter(wavenumber, isForward=False)
iopt=Optimiser(recTfunc, vdflt, isMonotonic=True, isAdaptive=True)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()

out=Output(flt, dflt, recTfunc, vdflt)
out.write()

k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', label="$G$")
plt.plot(k, dflt.getTransferFunction(), label="$G_d$")
plt.plot(k, recTfunc.getTransferFunction(), 'ro', label="${Q_d}^N$")
plt.plot(k, vgTfunc, label="$V_dG_d$")
plt.legend()
plt.grid('on')
plt.xlabel('kh')
plt.show()


