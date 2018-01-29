import numpy as np

def create_initial_function(Axis, mean, delta, sigma, momentum, dx=0.05):
    '''
    Return numpy array with the initial state of a wave function. The function is normalized.

    Parameters:
    Axis - array like
    mean - number, mean of gaussian shape wave
    delta - number, distance between two gaussians
    sigma - number, standard deviation of the gauss
    momentum - number, wave's momentum
    dx - float, optional, step in the axis
    '''

    InitialFunction = np.empty(len(Axis), dtype=complex)

    if delta == 0:
        for i, x in enumerate(Axis):
            InitialFunction[i] = np.exp(-(x-mean)*(x-mean)/(2*sigma*sigma))*np.exp(-momentum*x*1j)

    else:
        for i, x in enumerate(Axis):
            InitialFunction[i] = (np.exp(-(x-mean)*(x-mean)/(2*sigma*sigma))*np.exp(-momentum*x*1j)
                                  + np.exp(-(x-delta-mean)*(x-delta-mean)/(2*sigma*sigma))*np.exp(-momentum*x*1j))

                    
    #normalization          
    InitialFunction *= np.sqrt(1.0/(np.trapz(InitialFunction*InitialFunction.conjugate(), dx=dx))).real

    return InitialFunction


def create_absorption_filter(absorptionLength, dx=0.05):
    '''
    Return numpy array with absorption filter.

    Parameters:
    absorptionLength - int, width of the filter
    dx - float, optional, step in the axis
    '''

    AbsorptionFilter = np.empty(absorptionLength)

    for i in range(absorptionLength):
        if i > 0.1*absorptionLength:
            AbsorptionFilter[i] = np.exp(-(i*dx-int(absorptionLength*dx))*(i*dx-int(absorptionLength*dx))/3)
        else:
            AbsorptionFilter[i] = 0.0

    return AbsorptionFilter


def create_potential_function(potentialType, potentialWidth, Axis, dx=0.05):
    '''
    Return numpy array with potential function.

    Parameters:
    potentialType - char, shape of the barrier; 'r': rectangle,
                    't': isosceles triangle, 'rt1'/'rt2': right triangle 
    potentialWidth - number, width of the potential barrier
    Axis - array like    
    dx - float, optional, step in the axis
    '''

    axisLength = len(Axis)
    Potential = np.zeros(axisLength)
    middleX = (axisLength-1)/2

    #Rectangular potential barrier
    if potentialType == 'r':
        Potential[int(middleX-potentialWidth/(2*dx)):int(middleX+potentialWidth/(2*dx))] = 1


    #Triangular potential barrier
    if potentialType == 't':
        for i in range(int(middleX-potentialWidth/(2*dx)),int(middleX)):
            Potential[i] = (Axis[i] - Axis[int(middleX-potentialWidth/(2*dx))])/(Axis[int(middleX)]-Axis[int(middleX-potentialWidth/(2*dx))])

        for i in range(int(middleX),int(middleX+potentialWidth/(2*dx))):
            Potential[i] = (Axis[i] - Axis[int(middleX+potentialWidth/(2*dx))])/(Axis[int(middleX)]-Axis[int(middleX+potentialWidth/(2*dx))])


    if potentialType == 'rt1':
        for i in range(int(middleX-potentialWidth/(2*dx)),int(middleX+potentialWidth/(2*dx))):
            Potential[i] = (Axis[i] - Axis[int(middleX-potentialWidth/(2*dx))]) / potentialWidth


    if potentialType == 'rt2':
        for i in range(int(middleX-potentialWidth/(2*dx)),int(middleX+potentialWidth/(2*dx))):
            Potential[i] = (-Axis[i] - Axis[int(middleX-potentialWidth/(2*dx))]) / potentialWidth

        

    return Potential



def solve_hamiltonian(Potential, Axis, dx=0.05):
    '''
    For the given potential return eigenvalues and
    normalized eigenstates of a matrix representing Hamiltonian.

    Parameters:
    Potential - array like, potential distribution
    Axis - array like
    dx - float, step int the array
    '''

    axisLength = len(Axis)
    HamiltonianMatrix = np.empty(shape=(axisLength,axisLength))

    for i in range (axisLength):
        HamiltonianMatrix[i,i] = 1/(dx*dx) + Potential[i]

        try:
            HamiltonianMatrix[i,i+1] = -0.5/(dx*dx)
        except IndexError:
            pass

        try:
            HamiltonianMatrix[i,i-1] = -0.5/(dx*dx)
        except IndexError:
            pass
        

    Energy, StationaryStates = np.linalg.eigh(HamiltonianMatrix)

    # normalize stationary states
    for n in range(axisLength):
        StationaryStates[:,n] *= np.sqrt(1.0/(np.dot(StationaryStates[:,n],StationaryStates[:,n].conjugate()).real*dx))

    return Energy, StationaryStates



def create_coefficients(Function, StationaryStates, dx=0.05, numberOfStates=None):
    '''
    For the given function and stationary states
    return expansion coefficients and
    the sum of these coefficients.

    Parameters:
    Function - array like, wave function
    StationaryStates - array like, eigenvectors of the Hamiltionian
    dx - float, optional, step in the axis
    numberOfStates - int, optional, number of states in the expansion
    '''

    if numberOfStates is None:
        numberOfStates = StationaryStates.shape[1]

    Coefficients = np.empty(numberOfStates, dtype=complex)

    for i in range(numberOfStates):
        Coefficients[i] = np.trapz(Function*StationaryStates[:,i], dx=dx)
        
    coefficientsSum = np.dot(Coefficients, Coefficients.conjugate()).real

    return Coefficients, coefficientsSum


def create_function(t, axisLength, numberOfStates, Coefficients,
                   StationaryStates, Energy, TIME_COEFFICIENT):
    '''
    For the given time return the wave function.

    Parameters:
    t - number, time
    axisLength - int, length of the axis
    numberOfStates - int, number of states the expansion
    Coefficients - array like, coefficients from Fthe expansion
    SationaryStates - array like, eigenvalues of the Hamiltonian
    Energy - array like, energy of the stationary states
    TIME_COEFFICIENT - float, time coefficient
    '''

    return np.dot(StationaryStates[:,:numberOfStates],
                 (Coefficients[:numberOfStates]*np.exp(-Energy[:numberOfStates]*TIME_COEFFICIENT*t*1j)))


def absorption(Function, AbsorptionFilter, absorptionLength, axisLength, 
               transmissionCoefficient, reflectionCoefficient, dx=0.05):
    '''Return transmission and absorption coefficients.

    Parameters:
    Function - array like, wave function
    AbsorptionFilter - array like, absorption filter
    absorptionLength - int, width of the filter
    axisLength - inr, length of the axis
    transmissionCoefficient - float, transmission coefficient
    reflectionCoefficient - float, reflection coefficient
    dx - float, optional, step in the axis 
    '''

    transmission1 = np.trapz(Function[:absorptionLength]*Function[:absorptionLength].conjugate(), dx=dx).real
    reflection1 = np.trapz(Function[int(axisLength-absorptionLength):axisLength]*Function[int(axisLength-absorptionLength):axisLength].conjugate(), dx=dx).real

    Function[:absorptionLength] *= AbsorptionFilter[:absorptionLength]

    transmission2 = np.trapz(Function[:absorptionLength]*Function[:absorptionLength].conjugate(), dx=dx).real

    Function[(axisLength-absorptionLength):] *= AbsorptionFilter[absorptionLength::-1]

    reflection2 = np.trapz(Function[int(axisLength-absorptionLength):axisLength]*Function[int(axisLength-absorptionLength):axisLength].conjugate(), dx=dx).real

    transmissionCoefficient += transmission1-transmission2
    reflectionCoefficient += reflection1-reflection2

    return transmissionCoefficient, reflectionCoefficient

