from include import *
import scipy.signal as signal
import numpy as np

ks=np.pi*2.0
knyq=ks/2.0
deltak=0.10*knyq
kpass=np.pi/2.0-np.pi/8.0
kstop=np.pi/2.0+np.pi/8.0
ns=73
#
nsh=int((ns-1)/2)
#
wavenumber=Wavenumber()
k=wavenumber.get()
flt=Sharp(wavenumber, fType="bandpass", kbounds=[kpass, kstop])
df=DiscreteFilter(wavenumber)

#window method-Hamming
bw=signal.firwin(ns, [kpass, kstop], fs=ks, pass_zero=False)
hw=df.getTransferFunction(bw[nsh:])

#window method-Kaiser
bwk=signal.firwin(ns, [kpass, kstop], width=deltak, fs=ks, pass_zero=False)
hwk=df.getTransferFunction(bwk[nsh:])

#parks mcclellan
bp=signal.remez(ns, [0.0, kpass-deltak/2.0, kpass+deltak/2.0, kstop-deltak/2.0, kstop+deltak/2.0, np.pi], [0.0, 1.0, 0.0], fs=ks)
hp=df.getTransferFunction(bp[nsh:])

#least-squares optimisation
bo=signal.firls(ns, k, flt.getTransferFunction(), fs=ks)
fo, ho=signal.freqz(bo, fs=ks)
ho=df.getTransferFunction(bo[nsh:])

#pyDof
dflt=DiscreteFilter(wavenumber, nStencil=ns)
opt=Optimiser(flt, dflt, adaptive=False, fixZero=[0.0, np.pi], positive=True, lamda=30.0, error="mabse")  

opt.run()
out=Output(flt, dflt)
out.write()

#plt.plot(k, flt.getTransferFunction())
#plt.plot(k, flt.getTransferFunction(), 'k--', mfc="None")

plt.subplot(1,2,1)
#plt.plot(kc, ypass, 'k--')
#plt.plot(kstop, ypass, 'k--')
#plt.plot(k, (np.abs(hw)), label='$Window-Hamming$')
plt.plot(k, np.abs(hwk), label='$Window-Kaiser$')
plt.plot(k, np.abs(hp), label='$Parks-McClellan$')
plt.plot(k, np.abs(ho), label='$Least-squares$')
plt.plot(k, dflt.getTransferFunction(), label="$pyDof$")
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
#
plt.subplot(1,2,2)
#plt.plot(kc, 20*np.log10(ypass), 'k--')
#plt.plot(kstop, ypass, 'k--')
#plt.plot(k, (np.abs(hw)), label='$Window-Hamming$')
plt.plot(k, 20*np.log10(np.abs(hwk)), label='$Window-Kaiser$')
plt.plot(k, 20*np.log10(np.abs(hp)), label='$Parks-McClellan$')
plt.plot(k, 20*np.log10(np.abs(ho)), label='$Least-squares$')
plt.plot(k, 20*np.log10(dflt.getTransferFunction()), label="$pyDof$")
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-140, 0.05)
#
plt.show()



