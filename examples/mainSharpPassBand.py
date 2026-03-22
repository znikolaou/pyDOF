from include import *

nsList=[81]
wavenumber=Wavenumber()
k=wavenumber.get()
flt=SharpWindow(wavenumber, kcl=np.pi/4.0, kch=np.pi/2.0)
plt.plot(k, flt.getTransferFunction(), 'ko-', mfc="None", label="$\hat{G}$")
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False, isFullDampingAtK=[0.0, np.pi], isNoDampingAtK=[np.pi*3.0/8.0], lamda=3.0)
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

