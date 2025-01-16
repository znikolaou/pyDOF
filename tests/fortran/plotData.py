import pandas as pd
import matplotlib.pyplot as plt

data=pd.read_csv('data.txt')

plt.plot(data.get('y'),'b')
plt.plot(data.get('fy'),'r')
plt.plot(data.get('fya'),'ko')
plt.show()
