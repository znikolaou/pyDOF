from include import *

nsList=[81]
wavenumber=Wavenumber()
k=wavenumber.get()
flt=SharpWindow(wavenumber, kcl=np.pi/4.0, kch=np.pi/2.0)
plt.plot(k, flt.getTransferFunction(), 'ko-', mfc="None", label="$\hat{G}$")
maxPerc=[]
for v in nsList:
    m=int((v-1)/2)
    dflt=DiscreteFilter(wavenumber, nStencil=v)
    #opt=Optimiser(flt, dflt, isAdaptive=False)
    #opt=Optimiser(flt, dflt, isAdaptive=False, isFullDampingAtZero=True, isFullDampingAtPi=True, isPositive=True)
    opt=Optimiser(flt, dflt, isAdaptive=False, isZeroToOne=True, isFullDampingAtK=[0.0, np.pi], isNoDampingAtK=[np.pi*3.0/8.0])
    #opt=Optimiser(flt, dflt, isAdaptive=False, isFullDampingAtZero=True, isFullDampingAtPi=True, isZeroToOne=True, lamda=2.0)  
    opt.run()
    #mb=np.sqrt(np.dot(dflt.coeffs, dflt.coeffs))
    print("Coeffs=",dflt.coeffs)
    #print("cos(theta)=", np.dot(coeffsAnalytic, dflt.coeffs)/(ma*mb))
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

