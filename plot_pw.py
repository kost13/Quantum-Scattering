import numpy as np
from math import exp, e
from matplotlib import pyplot as plt
import os
        
#DATA FORM: 
#   0	  1			 2				 3			 4
# SIGMA ENERGY POTENTIAL_WIDTH  TRANSMISSION REFLECTION

ImportedData = np.loadtxt('../data/data_file.txt', delimiter=' ', dtype=float)


SIGMA = [r'$\sigma$', 0, 's']
ENERGY = [r'$ \cal E$', 1, 'e']
POTENTIAL_WIDTH = [r'$ V_0$' + ' width', 2, 'pw']
TRANSMISSION = ['Transmission', 3, 't']
REFLECTION = ['reflection', 4, 'r']

parameters = [SIGMA, ENERGY, POTENTIAL_WIDTH, TRANSMISSION, REFLECTION]

#SETTINGS

constant1 = ENERGY
value1 =  4
constant2 = SIGMA
Values = [0,0.5,2,10]

xData = POTENTIAL_WIDTH
yData = TRANSMISSION
functionAxis = plt.axes()

functionAxis.set_ylim([-0.005,1.005])
functionAxis.set_xlim([-0.005,8.005])

plt.grid(True)


directory = './graphs/'
dpi = 200


#END OF SETTINGS


col = xData[1]
ImportedData = ImportedData[np.argsort(ImportedData[:,col])]

functionAxis.set_xlabel(xData[0], fontsize=16)
functionAxis.set_ylabel(yData[0], fontsize=16)

style = ['bs', 'g^', 'ro', 'yp', 'cD', 'mv']

i = 0
for value in Values:
	CurveData = np.empty(shape=(0,5), dtype=float)
	for row in ImportedData:
		if row[constant1[1]] == value1 and row[constant2[1]] == value:
			CurveData = np.vstack([CurveData, row])

	X = CurveData[:,xData[1]]
	Y = CurveData[:,yData[1]]

	if value==0:
		plt.plot(X,Y, color='black', label='plain wave')
	else:
		plt.plot(X,Y, style[i], label=str(constant2[0]) + ' = ' + str(value), ms=4)
		i += 1


plt.legend(loc=4,numpoints=1)
if len(Values) == 1:
	FILE_NAME = './graphs/5-' + constant1[2] + str(round(value1,2)) + str(constant2[2]) + str(round(Values[0],2))
else:
	FILE_NAME = './graphs/5-' + constant1[2] + str(round(value1,2))


plt.tight_layout()

if os.path.isfile(FILE_NAME + '.png'):
    for i in range(1,1000):
        if not os.path.isfile(FILE_NAME + '(' + str(i) + ').png'):
            FILE_NAME = FILE_NAME + '(' + str(i) + ').png'
            break

else:
    FILE_NAME = FILE_NAME + '.png'

plt.savefig(FILE_NAME, dpi=dpi)
plt.show()