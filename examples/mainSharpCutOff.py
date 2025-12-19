from include import *

nsList=[11]

wavenumber=Wavenumber()
k=wavenumber.get()

flt=SharpCutOff(wavenumber, np.pi/4.0)
plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}$")
maxPerc=[]
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False)
    #opt=Optimiser(flt, dflt, isAdaptive=False, isNoDampingAtZero=True, isFullDampingAtPi=True, isZeroToOne=True)
    #opt=Optimiser(flt, dflt, isAdaptive=False, isNoDampingAtZero=True, isFullDampingAtPi=True, isMonotonicNegative=True)
    opt.run()
    print("Coeffs=",dflt.coeffs)
    maxPerc.append(getMaxPercDiff(flt, dflt))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, M=$"+str(m))
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.10, 1.1)
plt.show()
print('maxPerc=', maxPerc)

