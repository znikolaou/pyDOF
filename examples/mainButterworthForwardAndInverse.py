from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

wavenumber=Wavenumber()
flt=Butterworth(wavenumber, np.pi/2.0, 8.0)

dflt=DiscreteFilter(wavenumber)
opt=Optimiser(flt, dflt, optError=1.0e-4)
opt.run()

recTfunc=ReconstructedTransferFunction(dflt, nVanCittert=50)
vdflt=DiscreteFilter(wavenumber, isForward=False, nStencil=dflt.getNs())
iopt=Optimiser(recTfunc, vdflt, optError=1.0e-3)
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


