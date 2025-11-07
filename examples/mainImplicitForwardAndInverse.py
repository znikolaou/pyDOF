from include import *
from utils import *
from discreteFilter import *
from filters import *
from objectiveFunction import *
from optimiser import *
from output import * 

wavenumber=Wavenumber()
#flt=Implicit(wavenumber, -0.262)
flt=Implicit(wavenumber, -0.40)
dflt=DiscreteFilter(wavenumber)

opt=Optimiser(flt, dflt, isAdaptive=True)
opt.run()
ns=dflt.getNs()

recTfunc=ReconstructedTransferFunction(dflt, flt, nVanCittert=100)
vdflt=DiscreteFilter(wavenumber, nStencil=ns*2+1, isForward=False)
iopt=Optimiser(recTfunc, vdflt, isMonotonic=True, isAdaptive=False)
iopt.run()

vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()

out=Output("./output/", "coeffs", flt, dflt, recTfunc, vdflt)
out.write()

k=wavenumber.get()
plt.plot(k, flt.getTransferFunction(), 'ko', label="$G$")
plt.plot(k, dflt.getTransferFunction(), label="$G_d$")
plt.plot(k, recTfunc.getTransferFunction(), 'ro', label="${Q_d}^N$")
plt.plot(k, vgTfunc, label="$V_dG_d$")
plt.legend()
plt.grid('on')
plt.xlabel('$kh$')
plt.show()


