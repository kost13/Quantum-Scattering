#!/usr/bin/env python3
import sys
import warnings
import numpy as np
from matplotlib import pyplot as plt


def count_transmission_coefficient(a, E):
    '''
    Return coefficient of transmisson through potential barrier.
    Potential of the barrier equals 1.

    Parameters:
    a - width of the potential barrier
    E - energy of the wave
    '''

    if E == 1.0:
        E += 10e-7

    try:
        if E>1:
            return 1.0/(1 + (1 * (np.sin(np.sqrt(2*(E-1)) * a))**2)/(4*E*(E-1)))
        else:
            return 1.0/(1 + (1 * (np.sinh(np.sqrt(2*(1-E)) * a))**2)/(4*E*(1-E)))
    except RuntimeError:
        pass


def function_of_energy(potential_width = 4, energy_range = 5, step=0.01):
    '''
    Plot a function where y = Transmission, x = Energy of the wave

    parameters:
    potential_width - width of the barrier
    energy_range - range for x axis
    step - step for x axis
    '''

    Energy = np.arange(0,energy_range+step,step)
    TransmissionCoefficients = []

    for energy in Energy:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            functionValue = count_transmission_coefficient(potential_width, energy)

        TransmissionCoefficients.append(functionValue)

    plt.ylim([-0.01,1.01])
    plt.xlim([0,energy_range])
    plt.xlabel('Energy')
    plt.ylabel('Transmission')
    plt.grid(True)
    plt.plot(Energy, TransmissionCoefficients)
    plt.show()

def function_of_potential_width(energy = 1.1, potential_width_range = 8, step=0.01):
    '''
    Plot a function where y = Transmission, x = width of the potential barrier

    parameters:
    energy - energy of the wave
    potential_width_range - range of x axis
    step - step of x axis
    '''

    Potential_width = np.arange(0,potential_width_range+step,step)
    TransmissionCoefficients = []

    for potential_width in Potential_width:

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            functionValue = count_transmission_coefficient(potential_width, energy)

        TransmissionCoefficients.append(functionValue)

    plt.ylim([-0.01,1.01])
    plt.xlim([0,potential_width_range])
    plt.xlabel('Potential width')
    plt.ylabel('Transmission')
    plt.grid(True)
    plt.plot(Potential_width, TransmissionCoefficients)
    plt.show()

if __name__ == "__main__":

    function_of_energy()
    # function_of_potential_width()