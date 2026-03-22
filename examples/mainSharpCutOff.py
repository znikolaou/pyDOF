from include import *

nsList=[11]
coeffsAnalytic=[1/4, np.sqrt(2)/(2*np.pi), 1/(2*np.pi), np.sqrt(2)/(6*np.pi), 0.0, -np.sqrt(2)/(10*np.pi)]
ma=np.sqrt(np.dot(coeffsAnalytic, coeffsAnalytic))

wavenumber=Wavenumber()
k=wavenumber.get()

flt=SharpCutOff(wavenumber, np.pi/4.0)
plt.plot(k, flt.getTransferFunction(), 'ko-', mfc="None", label="$\hat{G}$")
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    opt=Optimiser(flt, dflt, isAdaptive=False)  
    opt.run()
    mb=np.sqrt(np.dot(dflt.coeffs, dflt.coeffs))
    print("Coeffs=",dflt.coeffs)
    print("cos(theta)=", np.dot(coeffsAnalytic, dflt.coeffs)/(ma*mb))
    out=Output(flt, dflt)
    out.write()
    plt.plot(k, dflt.getTransferFunction(), label="$\hat{G}_d, M=$"+str(m))
plt.xlabel('$kh$')
plt.legend(frameon=False, loc='center left')
plt.grid('on')
plt.xlim(0.0, np.pi)
plt.ylim(-0.05, 1.1)
plt.show()

