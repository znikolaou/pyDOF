from include import *

def getMaxOverShoot(dflt): return 100.0*(np.max(dflt.getTransferFunction())-1.0)

def getMaxPercDiff(flt, dflt, kc):
    return 100.0*np.max(np.abs(dflt.getTransferFunction()-flt.getTransferFunction())

nsList=[51, 81, 101]

wavenumber=Wavenumber()
k=wavenumber.get()

flt=SharpCutOff(wavenumber, np.pi/2.0)
plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}$")
maxVar=[]
for v in nsList:
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    #opt=Optimiser(flt, dflt, isAdaptive=False, isMonotonic=False, isZeroToOne=False, isFullDampingAtPi=False, isNoDampingAtZero=False, optError=v)
    opt=Optimiser(flt, dflt, isAdaptive=False, isMonotonic=True, isNoDampingAtZero=True, optError=v)
    opt.run()
    maxVar.append(getMaxOverShoot(dflt))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, ns=$"+str(v))
plt.xlabel('$kh$')
plt.legend(frameon=False)
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.05, 1.1)
plt.show()
print('maxVar=', maxVar)

