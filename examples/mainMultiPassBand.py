from include import *
import scipy.signal as signal
import numpy as np

ks=np.pi*2.0
knyq=ks/2.0
deltak=0.1*knyq
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

#Define a sharp passband filter
flt=Sharp(wavenumber, fType="bandpass", kbounds=bands)
df=DiscreteFilter(wavenumber)
tfSharp=flt.getTransferFunction()

#window method-Kaiser
bwk=signal.firwin(ns, bands, width=deltak, window='kaiser', fs=ks, pass_zero='bandpass')
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
r=0.0005
kStopA=[0, ka-deltak/2.0]
kPassA=[ka+deltak/2.0, kb-deltak/2.0]
kStopB=[kb+deltak/2.0, kc-deltak/2.0]
kPassB=[kc+deltak/2.0, kd-deltak/2.0]
kStopC=[kd+deltak/2.0, np.pi]
opt=Optimiser(flt, dflt, adaptive=False, positive=False, fixZero=[0.0, np.pi], \
        luBounds=[ [-r, r, kStopA], [-r, r, kStopB], [-r, r, kStopC], [1-r, 1+r, kPassA], [1-r, 1+r, kPassB] ], error="msqe")  
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

plt.plot(k, tfSharp, label='$Ideal$')
plt.plot(k, (np.abs(hwk)), label='$Window-Kaiser$')
plt.plot(k, (np.abs(hp)), label='$Parks-McClellan$')
plt.plot(k, (np.abs(ho)), label='$Least-squares$')
plt.plot(k, (dflt.getTransferFunction()), '--', label="$pyDof$")
plt.xlabel('$kh$')
plt.ylabel('$\hat{G}_d$')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.1, 1.1)
plt.show()

plt.plot(k, 20*np.log10(np.abs(hwk)), label='$Window-Kaiser$')
plt.plot(k, 20*np.log10(np.abs(hp)), label='$Parks-McClellan$')
plt.plot(k, 20*np.log10(np.abs(ho)), label='$Least-squares$')
plt.plot(k, 20*np.log10(dflt.getTransferFunction()), '--', label="$pyDof$")
plt.xlabel('$kh$')
plt.ylabel('$dB(\hat{G}_d)$')
plt.legend(frameon=False, loc='lower left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-130, 10.0)
plt.show()



