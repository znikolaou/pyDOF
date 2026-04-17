from include import *
from errors import *
import scipy.signal as signal
import numpy as np

ks=np.pi*2.0
knyq=ks/2.0
deltak=0.10*knyq
kcut=knyq/2.0
ns=73
#
nsh=int((ns-1)/2)
#
wavenumber=Wavenumber()
k=wavenumber.get()
ipass=k.searchsorted(kcut-deltak/2.0)
kpass=k[0:ipass+1]
df=DiscreteFilter(wavenumber)

#Define a sharp lowpass filter
flt=Sharp(wavenumber, kbounds=[kcut], fType='lowpass')
tfSharp=flt.getTransferFunction()

#window method-Hamming
bw=signal.firwin(ns, cutoff=kcut, fs=ks)
hw=df.getTransferFunction(bw[nsh:])

#window method-Kaiser
bwk=signal.firwin(ns, cutoff=kcut, width=deltak, window='kaiser', fs=ks, pass_zero='lowpass')
hwk=df.getTransferFunction(bwk[nsh:])

#parks mcclellan
bp=signal.remez(ns, [0.0, kcut-deltak/2.0, kcut+deltak/2.0, np.pi], [1.0, 0.0], fs=ks)
hp=df.getTransferFunction(bp[nsh:])

#least-squares optimisation
bo=signal.firls(ns, k, flt.getTransferFunction(), fs=ks)
fo, ho=signal.freqz(bo, fs=ks)
ho=df.getTransferFunction(bo[nsh:])

#pyDof
dflt=DiscreteFilter(wavenumber, nStencil=ns)
#opt=Optimiser(flt, dflt, isAdaptive=False, isPositive=False, lamda=20.0, error="mabse")  
opt=Optimiser(flt, dflt, isAdaptive=False, isPositive=False, isNoDampingAtK=[0.0, kcut-deltak/2.0], isFullDampingAtK=[kcut+deltak/2.0, np.pi], isBounded=[[1-0.0005, 1+0.0005, [0, kcut-deltak/2.0]], [-0.0005, 0.0005, [kcut+deltak/2.0, np.pi]]], lamda=0.0, error="msqe")
#opt=Optimiser(flt, dflt, isAdaptive=False, isPositive=False, isBounded=[1-0.0005, 1+0.0005, kpass], lamda=10.0, error="mabse")  
#
opt.run()
#out=Output(flt, dflt)
#out.write()

#errors
msqe={}
msqe['WK']=Errors.msqe(hwk, tfSharp)
msqe['PM']=Errors.msqe(hp, tfSharp)
msqe['SQ']=Errors.msqe(ho, tfSharp)
msqe['PDA']=Errors.msqe(dflt.getTransferFunction(), tfSharp)

mabse={}
mabse['WK']=Errors.mabse(hwk, tfSharp)
mabse['PM']=Errors.mabse(hp, tfSharp)
mabse['SQ']=Errors.mabse(ho, tfSharp)
mabse['PDA']=Errors.mabse(dflt.getTransferFunction(), tfSharp)

#max overshoot
overshoot={}
overshoot['WK']=(np.max(hwk)-1.0)*100.0
overshoot['PM']=(np.max(hp)-1.0)*100.0
overshoot['SQ']=(np.max(ho)-1.0)*100.0
overshoot['PDA']=(np.max(dflt.getTransferFunction())-1.0)*100.0

print('MSQE:')
for key in msqe: print(key, msqe[key])
print('MABSE:')
for key in mabse: print(key, mabse[key])
print('Max % overshoot in passband:')
for key in overshoot: print(key, overshoot[key])

plt.subplot(1,2,1)
plt.plot(k, tfSharp, 'k--')
plt.plot(k, hwk, label='$Window-Kaiser$')
plt.plot(k, hp, label='$Parks-McClellan$')
plt.plot(k, ho, label='$Least-squares$')
plt.plot(k, dflt.getTransferFunction(), '--', label="$pyDof-A$")
plt.xlabel('$kh$')
plt.ylabel('$\hat{G}_d$')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)

plt.subplot(1,2,2)
plt.plot(k, 20*np.log10(hwk), label='$Window-Kaiser$')
plt.plot(k, 20*np.log10(hp), label='$Parks-McClellan$')
plt.plot(k, 20*np.log10(ho), label='$Least-squares$')
plt.plot(k, 20*np.log10(dflt.getTransferFunction()), '--', label="$pyDof-A$")
plt.xlabel('$kh$')
plt.ylabel('dB')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-140, 5.0)
#
plt.show()



