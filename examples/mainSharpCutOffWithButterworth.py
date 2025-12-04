from include import *

def getMaxOverShoot(dflt): return 100.0*(np.max(dflt.getTransferFunction())-1.0)

nsList=[81]

wavenumber=Wavenumber()
k=wavenumber.get()

flt=SharpCutOff(wavenumber, np.pi/4.0)
fltB=Butterworth(wavenumber, np.pi/4.0, 6)

plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}$")
plt.plot(k, fltB.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}$")
maxVar=[]
maxVarB=[]
for v in nsList:
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    dfltB=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False, isMonotonic=False, isPositive=False, isZeroDampingAtPi=False, optError=v)
    optB=Optimiser(fltB, dfltB, isAdaptive=False, isMonotonic=False, isPositive=False, isZeroDampingAtPi=False, optError=v)
    opt.run()
    optB.run()
    maxVar.append(getMaxOverShoot(dflt))
    maxVarB.append(getMaxOverShoot(dfltB))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, ns=$"+str(v))
    plt.plot(k, dfltB.getTransferFunction(), label="$\hat{G}_d-B, ns=$"+str(v))
plt.xlabel('$kh$')
plt.legend(frameon=False)
plt.grid('on')
plt.show()
print('maxVar=', maxVar)
print('maxVarB=', maxVarB)

