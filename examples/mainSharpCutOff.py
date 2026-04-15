from include import *
import scipy.signal as signal
import numpy as np

fs=np.pi*2.0
piOv2=np.pi/2.0
ns=51
transWidth=0.10*piOv2
#
nsh=int((ns-1)/2)

wavenumber=Wavenumber()
k=wavenumber.get()
flt=SharpCutOff(wavenumber, piOv2)
df=DiscreteFilter(wavenumber)

#window method
bw=signal.firwin(ns, cutoff=piOv2, fs=fs)
#fw, hw=signal.freqz(bw, fs=fs)
hw=df.getTransferFunction(bw[nsh:])

#parks mclelan
bp=signal.remez(ns, [0.0, piOv2-transWidth/2.0, piOv2+transWidth/2.0, np.pi], [1.0, 0.0], fs=fs)
#fp, hp=signal.freqz(bp, fs=fs)
hp=df.getTransferFunction(bp[nsh:])

#least-squares optimisation
bo=signal.firls(ns, k, flt.getTransferFunction(), fs=fs)
fo, ho=signal.freqz(bo, fs=fs)
ho=df.getTransferFunction(bo[nsh:])

nsList=[ns]
#coeffsAnalytic=[1/4, np.sqrt(2)/(2*np.pi), 1/(2*np.pi), np.sqrt(2)/(6*np.pi), 0.0, -np.sqrt(2)/(10*np.pi)]
#ma=np.sqrt(np.dot(coeffsAnalytic, coeffsAnalytic))

plt.plot(k, flt.getTransferFunction(), 'ko-', mfc="None", label="$\hat{G}$")
plt.plot(k, np.abs(hw), label='$Hamming window$')
plt.plot(k, np.abs(hp), label='$Parks-McLelan$')
plt.plot(k, np.abs(ho), label='$LS optimiser: firls$')
plt.legend()
#plt.show()
#exit()
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False, isNoDampingAtK=[0.0], isPositive=True, lamda=3.0)  
    #opt=Optimiser(flt, dflt, isAdaptive=False)  
    
    opt.run()
    #mb=np.sqrt(np.dot(dflt.coeffs, dflt.coeffs))
    #print("Coeffs=",dflt.coeffs)
    #print("cos(theta)=", np.dot(coeffsAnalytic, dflt.coeffs)/(ma*mb))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$pyDof:\hat{G}_d, M=$"+str(m))

plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.05, 1.1)
plt.show()

c=opt.get()
for i in range(nsh):
    print(i, bw[i+nsh], bp[i+nsh], bo[i+nsh], c[i])


