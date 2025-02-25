import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

def readBinary(nx, ny, file):
    integ='i4'
    prec='f8'
    f=FortranFile(file)
    data=f.read_record( \
            integ, \
            integ, \
            np.dtype((prec,(ny,nx))), \
            np.dtype((prec,(ny,nx))), \
            np.dtype((prec,(ny,nx))) ) 
    nx=data[0]
    ny=data[1]
    a=data[2]
    fa=data[3]
    ra=data[4]
    return a, fa, ra

a, fa, ra=readBinary(400, 400, 'data.bin')

plt.subplot(4,1,1)
plt.contourf(a, label='a')
plt.colorbar()
plt.subplot(4,1,2)
plt.contourf(fa, label='fa')
plt.colorbar()
plt.subplot(4,1,3)
plt.contourf(ra, label='ra')
plt.colorbar()
plt.subplot(4,1,4)
plt.contourf(ra-a, label='ra-a')
plt.colorbar()
plt.show()
