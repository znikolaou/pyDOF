import pandas as pd
import matplotlib.pyplot as plt

data=pd.read_csv('data.txt')

plt.plot(data.get('y'),'b', label='y')
plt.plot(data.get('fy'),'r', label='fy')
plt.plot(data.get('fya'),'ko', label='fy-actual')
plt.plot(data.get('ry'),'go', label='ry')
plt.legend()
plt.show()
