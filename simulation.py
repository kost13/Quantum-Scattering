'''
Scattering of wave function.
'''

import numpy as np
from functions import *

#EDIT

SIGMA = [2]

ENERGY = np.arange(0.1, 4.6, 0.1)
POTENTIAL_BARRIER_WIDTH = [1,2,4,8]

#axis settings
step = 0.05
xRange = 25

mean = 8
delta = 0


maxTimeLength = 80
timeCoefficient = 0.05

potentialType = 'r'

#default number of states, this value if self-adjusting
numberOfStates=200

# END OF EDIT


absorptionLength = int(5/step)

Axis = np.arange(-xRange-absorptionLength*step, xRange+absorptionLength*step+step, step)
axisLength=len(Axis)   


AbsorptionFilter = create_absorption_filter(absorptionLength, step)


def main(sigma, potentialWidth, energy, FILE_NAME):

	momentum = np.sqrt(2.0*energy)

	InitialFunction = create_initial_function(Axis, mean, delta, sigma, momentum, step)

	Potential = create_potential_function(potentialType, potentialWidth, Axis, step)
	
	Energy, StationaryStates = solve_hamiltonian(Potential, Axis, step)

	reflectionCoefficient = 0.0
	transmissionCoefficient = 0.0

	Function = InitialFunction
	
	Coefficients, coefficentsSum  = create_coefficients(Function, StationaryStates, step, numberOfStates)

	coefficientsAccuracy = coefficentsSum + reflectionCoefficient + transmissionCoefficient

	coefficentsRepeat = 0
	time = 0

	# conditions to stop symulation 
	while(round(transmissionCoefficient + reflectionCoefficient,3) < 0.999
		  and round(transmissionCoefficient,3) + round(reflectionCoefficient,3) < 0.999
		  and time < maxTimeLength/timeCoefficient
		  and (coefficentsRepeat < 8 or (transmissionCoefficient + reflectionCoefficient) < 0.99)
		  ):

		Function = create_function(time, axisLength, numberOfStates, Coefficients,
			                       StationaryStates, Energy, timeCoefficient)

		transmissionCoefficient, reflectionCoefficient = absorption(Function, AbsorptionFilter, absorptionLength,axisLength,
			                                                        transmissionCoefficient, reflectionCoefficient, step)

		lastCoefficientsAccuracy = coefficientsAccuracy

		Coefficients, coefficentsSum = create_coefficients(Function, StationaryStates, step, numberOfStates)
		coefficientsAccuracy = coefficentsSum + reflectionCoefficient + transmissionCoefficient

		#checks if further symulation can improve the result
		if abs(coefficientsAccuracy-lastCoefficientsAccuracy) < 10e-9:
			coefficentsRepeat += 1
		else:
			coefficentsRepeat = 0

		if coefficientsAccuracy < 0.9991:
			return False

		if (time % 10) == 0:
			print(str(timeCoefficient*time) + '  ' + str(numberOfStates) + ' ' + str(coefficientsAccuracy) + ' '
				  + str(round(transmissionCoefficient,3)) +  '  ' + str(round(reflectionCoefficient,3)) + '  '
				  + str(round(transmissionCoefficient+reflectionCoefficient,3)))
		time += 1


	print(str(timeCoefficient*time)  + '  ' + str(numberOfStates) + ' '
		  + str(coefficientsAccuracy) + '  ' + str(round(transmissionCoefficient,3)) +  '  '
		  + str(round(reflectionCoefficient,3)) + '  '  + str(round(transmissionCoefficient+reflectionCoefficient,3)))
	
	print('\n')

	data = (str(sigma) + ' ' + str(energy) + ' '
		    + str(round(potentialWidth,4)) + ' ' + str(round(transmissionCoefficient,3)) + ' '
		    + str(round(reflectionCoefficient,3)) + '\n')

	# data is saved in .txt file in folder data
	with open('data/' + str(FILE_NAME), 'a') as file:
		file.write(data)

	return True

for sigma in SIGMA:
	for potentialWidth in POTENTIAL_BARRIER_WIDTH:
		for energy in ENERGY:

			FILE_NAME = "simulation_data.txt"

			print "E: ", energy, " pw: ", potentialWidth, " sigma: ", sigma

			# if it is impossible to get high precision with current
			# number of states, the number of states increases.
			ok = False
			while not ok:
				ok = main(sigma, potentialWidth, energy, FILE_NAME)

				if not ok:
					print "repeat"
					numberOfStates = min(int(numberOfStates*1.2), axisLength)			