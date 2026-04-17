from include import *
import scipy.signal as signal
import numpy as np

ks=np.pi*2.0
knyq=ks/2.0
deltak=0.10*knyq
ka=np.pi/4.0-np.pi/8.0/2.0
kb=np.pi/4.0+np.pi/8.0/2.0
kc=np.pi*3.0/4.0-np.pi/8.0/2.0
kd=np.pi*3.0/4.0+np.pi/8.0/2.0
ns=73
#
bands=[ka, kb, kc, kd]
nsh=int((ns-1)/2)
#
wavenumber=Wavenumber()
k=wavenumber.get()
flt=Sharp(wavenumber, fType="bandpass", kbounds=bands)
df=DiscreteFilter(wavenumber)

#window method-Hamming
#bw=signal.firwin(ns, bands, fs=ks, pass_zero=False)
#hw=df.getTransferFunction(bw[nsh:])

#window method-Kaiser
bwk=signal.firwin(ns, bands, width=deltak, fs=ks, pass_zero=False)
hwk=df.getTransferFunction(bwk[nsh:])

#parks mcclellan
bp=signal.remez(ns, [0.0, ka-deltak/2.0, ka+deltak/2.0, kb-deltak/2.0, kb+deltak/2.0, kc-deltak/2.0, kc+deltak/2.0, kd-deltak/2.0, kd+deltak/2.0, np.pi], [0.0, 1.0, 0.0, 1.0, 0.0], fs=ks)
hp=df.getTransferFunction(bp[nsh:])

#least-squares optimisation
bo=signal.firls(ns, k, flt.getTransferFunction(), fs=ks)
fo, ho=signal.freqz(bo, fs=ks)
ho=df.getTransferFunction(bo[nsh:])

#pyDof
dflt=DiscreteFilter(wavenumber, nStencil=ns)

#works very well!
opt=Optimiser(flt, dflt, isAdaptive=False, isFullDampingAtK=[0.0, np.pi], isNoDampingAtK=[np.pi/4.0, np.pi*3.0/4.0], isPositive=True,  lamda=30.0, error="mabse")  
opt.run()
#out=Output(flt, dflt)
#out.write()

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



