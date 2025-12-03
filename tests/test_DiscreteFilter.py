from include import *
from discreteFilter import *
from utils import * 
import pytest

waveno=Wavenumber()

def testIsSymmetric():
    with pytest.raises(SystemExit) as e:
        dflt=DiscreteFilter(waveno, nStencil=8)
    assert e.type==SystemExit
  
def testIsSufficientStencil():
    with pytest.raises(SystemExit) as e:
        dflt=DiscreteFilter(waveno, nStencil=1)
    assert e.type==SystemExit
  


