from scipy import signal
import matplotlib.pyplot as plt 
import numpy as np

sos=signal.butter(N=4, Wn=2.0, btype='lowpass', fs=10, output='sos')
b, a=signal.butter(N=4, Wn=2.0, btype='lowpass', fs=10, output='sos')

print(b, a)

w,h=signal.sosfreqz(sos, worN=8000)

plt.plot(w,np.abs(h), 'ko')
plt.show()

