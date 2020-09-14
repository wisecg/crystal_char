import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

table = pd.read_csv('/Users/ccenpa/mytemp.txt', delimiter=' ')
index = np.arange(0,628)
near = table['Near']
far = table['Far']

# print(table['Time'][191])

def lin(x, *p):
    a,b=p
    return a*x + b

p1=[1,1]
popt1,pcov1=curve_fit(lin, index[191:], near[191:], p0=p1)
print(popt1)

p2=[1,1]
popt2,pcov2=curve_fit(lin, index[191:], far[191:], p0=p2)
print(popt2)

fit1 = lin(index[191:], *popt1)
fit2 = lin(index[191:], *popt2)

# plt.plot(index, near, c='xkcd:aqua', marker='.', linestyle='none', label='Temperature close to HV')
# plt.plot(index, far, c='xkcd:yellow orange', marker='.', linestyle='none', label='Temperature far from HV')
# plt.plot(index[191:], fit1, c='xkcd:aqua')
# plt.plot(index[191:], fit2, c='xkcd:yellow orange')
# plt.xlabel('Time (~10 second intervals)')
# plt.ylabel('Temperature, degrees Celsius')
# plt.legend()
# plt.show()

rvx = pd.read_csv('/Users/ccenpa/roomvxtal.txt', delimiter=' ')

index2 = np.arange(0, 1427)

plt.plot(index2, rvx['In'], c='xkcd:sky', marker='.', linestyle='none', label='Temperature inside crystal')
plt.plot(index2, rvx['Out'], c='xkcd:deep red', marker='.', linestyle='none', label='Temperature outside crystal')
# plt.plot(index[191:], fit1, c='xkcd:aqua')
# plt.plot(index[191:], fit2, c='xkcd:yellow orange')
plt.xlabel('Time (~10 second intervals)')
plt.ylabel('Temperature, degrees Celsius')
plt.legend()
plt.show()
