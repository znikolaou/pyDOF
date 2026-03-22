import matplotlib.pyplot as plt
import numpy as np
#----------------------------------------------------------------------
small_size=24
medium_size=24
big_size=24
legend_size=20
plt.rc('font', size=small_size)         
plt.rc('axes', titlesize=big_size)      
plt.rc('axes', labelsize=big_size)       
plt.rc('xtick', labelsize=big_size)     
plt.rc('ytick', labelsize=big_size)      
plt.rc('legend', fontsize=legend_size)  
plt.rc('figure', titlesize=big_size)  
plt.rcParams['figure.figsize']=(11,8)
plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

def toDb(x): return 20.0*np.log10(x)

def plotAll(k, flt, dflt, recTfunc, vdflt, nVanCit):
    vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
    plt.plot(k, flt.getTransferFunction(), 'ko', mfc="None", markevery=5, label="$\hat{G}$")
    plt.plot(k, dflt.getTransferFunction(), 'r', label="$\hat{G}_d$")
    plt.plot(k, recTfunc.getTransferFunction(), 'bd', mfc="None", markevery=5, label="${\hat{Q}_d}^N$")
    plt.plot(k, vdflt.getTransferFunction()/np.max(vdflt.getTransferFunction()), 'g', label="$\hat{V}_d+ (N=$"+str(nVanCit)+")")
    plt.plot(k, vgTfunc, 'm', label="$\hat{V}_d\hat{G}_d$")
    plt.xlim(0, np.pi)
    plt.ylim(-0.05, 1.1)
    plt.xlabel('$kh$')
    plt.legend(frameon=False, loc='center left')
    plt.grid('on')
    plt.show()

def plotAllDb(k, flt, dflt, recTfunc, vdflt, nVanCit):
    vgTfunc=vdflt.getTransferFunction()*dflt.getTransferFunction()
    plt.plot(k, toDb(flt.getTransferFunction()), 'ko', mfc="None", markevery=5, label="$\hat{G}$")
    plt.plot(k, toDb(dflt.getTransferFunction()), 'r', label="$\hat{G}_d$")
    plt.plot(k, toDb(recTfunc.getTransferFunction()), 'bd', mfc="None", markevery=5, label="${\hat{Q}_d}^N$")
    plt.plot(k, toDb(vdflt.getTransferFunction()/np.max(vdflt.getTransferFunction())), 'g', label="$\hat{V}_d+ (N=$"+str(nVanCit)+")")
    plt.plot(k, toDb(vgTfunc), 'm', label="$\hat{V}_d\hat{G}_d$")
    plt.xlim(0, np.pi)
    plt.ylim(-60.0, 10.0)
    plt.xlabel('$kh$')
    #plt.ylabel('$dB$')
    plt.legend(frameon=False, loc='center left')
    plt.grid('on')
    plt.show()

