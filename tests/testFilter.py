from include import *
from filters import *
import numpy as np

k=Wavenumber()
butter=Butterworth(k, 4, np.pi/4.0)

plt.plot(k.get(), butter.getTransferFunction())
plt.show()


