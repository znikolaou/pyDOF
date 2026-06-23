from include import *

ks=np.pi*2.0
knyq=ks/2.0
deltak=0.05*knyq
kcut=0.5*knyq
ns=145
rpl=0.0005

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
opt=Optimiser(flt, dflt, adaptive=False, fixOne=[0.0], fixZero=[np.pi], luBounds=[[1-rpl, 1+rpl, [0, kcut-deltak/2.0]], [-rpl, rpl, [kcut+deltak/2.0, np.pi]]])
opt.run()
out=Output(flt, dflt)
out.write()

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

plt.plot(k, (hwk), label='$Kaiser$')
plt.plot(k, (hp), label='$Parks-McClellan$')
plt.plot(k, (ho), label='$Least-squares$')
plt.plot(k, (dflt.getTransferFunction()), '--', label="$pyDof$")
plt.xlabel('$kh$')
plt.ylabel('$\hat{G}_d$')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.1, 1.1)
plt.show()

plt.plot(k, 20*np.log10(hwk), label='$Kaiser$')
plt.plot(k, 20*np.log10(hp), label='$Parks-McClellan$')
plt.plot(k, 20*np.log10(ho), label='$Least-squares$')
plt.plot(k, 20*np.log10(dflt.getTransferFunction()), '--', label="$pyDof$")
plt.xlabel('$kh$')
plt.ylabel('$dB(\hat{G}_d)$')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-130, 10.0)
plt.show()


