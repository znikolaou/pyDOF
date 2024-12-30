from include import *
from filters import *

k=Wavenumber()
gaus=Gaussian(k, 4.0)
plt.plot(k.get(), gaus.get())
plt.show()

print(gaus.getCutOff())
print(gaus.getNyqDamping())
print(gaus.getEfficiency())


