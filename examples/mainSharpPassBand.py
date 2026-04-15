from include import *
import scipy.signal as signal

fs=np.pi*2.0
nsList=[51]
nsh=int((nsList[0]-1)/2)

wavenumber=Wavenumber()
k=wavenumber.get()
flt=SharpWindow(wavenumber, kcl=np.pi/4.0, kch=np.pi/2.0)
df=DiscreteFilter(wavenumber)

#window method
bw=signal.firwin(nsList[0], [np.pi/4.0, np.pi/2.0], fs=fs, pass_zero=False)
#fw, hw=signal.freqz(bw, fs=np.pi*2.0)
hw=df.getTransferFunction(bw[nsh:])

#parks mclelan
transWidth=0.1*np.pi/2.0
bp=signal.remez(nsList[0], [0.0, np.pi/4.0-transWidth/2.0, np.pi/4.0+transWidth/2.0, np.pi/2.0-transWidth/2.0, np.pi/2.0+transWidth/2.0, np.pi], [0.0, 1.0, 0.0], fs=fs)
#fp, hp=signal.freqz(bp, fs=fs)
hp=df.getTransferFunction(bp[nsh:])

#least squares
bo=signal.firls(nsList[0], k, flt.getTransferFunction(), fs=fs)
#fo, ho=signal.freqz(bp, fs=np.pi*2.0)
ho=df.getTransferFunction(bo[nsh:])

plt.plot(k, flt.getTransferFunction(), 'ko-', mfc="None", label="$\hat{G}$")
plt.plot(k, hw, label="window")
plt.plot(k, hp, label="Parks Mclelan")
plt.plot(k, ho, label="firls")
#plt.show()
#exit()
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False, isFullDampingAtK=[0, np.pi], isNoDampingAtK=[np.pi/4.0+transWidth/2.0, np.pi*3.0/8.0, np.pi/2.0-transWidth/2.0], lamda=5.5)
    opt.run()
    print("Coeffs=",dflt.coeffs)
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, M=$"+str(m))
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.05, 1.1)
plt.show()

