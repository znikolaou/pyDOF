from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

wavenumber=Wavenumber()
flt=Butterworth(wavenumber, np.pi/2.0, 6.0)

dflt=DiscreteFilter(wavenumber, nStencil=31)
opt=Optimiser(flt, dflt, isZeroDampingAtPi=False, isMonotonic=True, isAdaptive=False)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=50)
vdflt=DiscreteFilter(wavenumber, isForward=False, nStencil=31)
iopt=Optimiser(recTfunc, vdflt, isMonotonic=True, isAdaptive=False)
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


