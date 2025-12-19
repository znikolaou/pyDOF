from include import *

def getMaxPercDiff(flt, dflt):
    iPass=flt.getPassBandIndex()
    pPass=100.0*(np.max(dflt.getTransferFunction()[0:iPass+1]-1.0))
    return pPass

nsList=[161]
nB=32

wavenumber=Wavenumber()
k=wavenumber.get()

flt=SharpCutOff(wavenumber, np.pi/2.0)
fltB=Butterworth(wavenumber, np.pi/2.0, nB)
plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", label="$\hat{G}-sharp cut-off$")
plt.plot(k, fltB.getTransferFunction(), 'k--', mfc="None", label="$\hat{G}-Butterworth (n="+str(nB)+")$")

maxPerc=[]
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(fltB, dflt, isAdaptive=False, isMonotonic=True, isZeroToOne=False,  isNoDampingAtZero=True, isFullDampingAtPi=True, optError=v)
    opt.run()
    maxPerc.append(getMaxPercDiff(flt, dflt))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, M=$"+str(m))
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.05, 1.1)
plt.show()
print('maxPerc=', maxPerc)

