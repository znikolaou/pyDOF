from include import *
from utils import *

waveno=Wavenumber()

def testK():
    assert waveno.get().all()==np.linspace(0.0e0, np.pi, waveno.getNk()).all()
