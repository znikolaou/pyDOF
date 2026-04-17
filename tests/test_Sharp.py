from include import *
from utils import * 
from filters import *
import pytest

waveno=Wavenumber(20)

def testExitSizeLowPass():
    with pytest.raises(SystemExit) as e:
        flt=Sharp(waveno, fType="lowpass", kbounds=[1, 2])
    assert e.type==SystemExit

def testExitSizeBandPass():
    with pytest.raises(SystemExit) as e:
        flt=Sharp(waveno, fType="bandpass", kbounds=[1, 2, 3])
    assert e.type==SystemExit

def testExitBoundsHigh():
    with pytest.raises(SystemExit) as e:
        flt=Sharp(waveno, fType="lowpass", kbounds=[2.0*np.pi])
    assert e.type==SystemExit

def testExitBoundsLow():
    with pytest.raises(SystemExit) as e:
        flt=Sharp(waveno, fType="lowpass", kbounds=[-2.0*np.pi])
    assert e.type==SystemExit

def testExitNonMonotone():
    with pytest.raises(SystemExit) as e:
        flt=Sharp(waveno, fType="bandpass", kbounds=[2, 1])
    assert e.type==SystemExit

def testLowPass():
    flt=Sharp(waveno, fType="lowpass", kbounds=[np.pi/2.0])
    k=waveno.get()
    tf=flt.getTransferFunction()
    for i in range(len(k)): print(k[i], tf[i])

def testHighPass():
    flt=Sharp(waveno, fType="highpass", kbounds=[np.pi/2.0])
    k=waveno.get()
    tf=flt.getTransferFunction()
    for i in range(len(k)): print(k[i], tf[i])

def testBandPass():
    flt=Sharp(waveno, fType="bandpass", kbounds=[np.pi/4.0, np.pi/2.0])
    k=waveno.get()
    tf=flt.getTransferFunction()
    for i in range(len(k)): print(k[i], tf[i])

def testBandStop():
    flt=Sharp(waveno, fType="bandstop", kbounds=[np.pi/4.0, np.pi/2.0])
    k=waveno.get()
    tf=flt.getTransferFunction()
    for i in range(len(k)): print(k[i], tf[i])

def testMultiBandPass():
    flt=Sharp(waveno, fType="bandpass", kbounds=[np.pi/8.0, np.pi/4.0, np.pi*3.0/8.0, np.pi/2.0, np.pi*5.0/8.0, np.pi*6.0/8.0])
    k=waveno.get()
    tf=flt.getTransferFunction()
    for i in range(len(k)): print(k[i], tf[i])




  

