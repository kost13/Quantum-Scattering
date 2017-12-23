import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import matplotlib.lines as mlines
import os

from functions import *

sigma = 2
energy = 2
potentialWidth = 4


#axis settings
step = 0.05
xFrom = -20
xTo = 20

mean = 8
delta = 0

timeCoefficient = 0.1

potentialType = 'rt2'

yRange = 0.4

potentialDisplay = True 
gridDisplay = True 
#minimum value: 0.000001
timeConverter = 0.01

#time of disaplying one frame in ms
intervalTime = 100

fpsNumber = 10
animationName = 'aaaa-gauss6'
directory = './animations/'
bitRate = 1200

# ANIMATION_DURATION * intervalTime = animation length in ms
ANIMATION_DURATION = 120



absorptionLength = int(5/step)

xlim = (xFrom+5,xTo-5)

Axis = np.arange(xFrom, xTo+step, step)
axisLength=len(Axis)   

numberOfStates=400

momentum = np.sqrt(2*energy)

AbsorptionFilter = create_absorption_filter(absorptionLength, step)


InitialFunction = create_initial_function(Axis, mean, delta, sigma, momentum, step)

Potential = create_potential_function(potentialType, potentialWidth, Axis, step)

Energy, StationaryStates = solve_hamiltonian(Potential, Axis, step)

reflectionCoefficient = 0.0
transmissionCoefficient = 0.0

Function = InitialFunction

Coefficients, coefficentsSum  = create_coefficients(Function, StationaryStates, step, numberOfStates)



fig = plt.figure()
functionAxis = plt.axes(xlim=xlim, ylim=(0,yRange))

functionLine, = functionAxis.plot([], [], lw=1)

functionAxis.set_xlabel('x', fontsize=16)
functionAxis.set_ylabel(r'$|\psi|^2$', fontsize=18, rotation='horizontal', ha='right')

plt.tight_layout()
fig.subplots_adjust(top=.90)

props = dict(boxstyle='square, pad=0.5', facecolor='white')
timeText = functionAxis.text(0.012, 1.05, '', transform=functionAxis.transAxes, fontsize=13, bbox=props)


def init():
    functionLine.set_data([], [])
    timeText.set_text('')
    return functionLine, timeText


def animate(i):
    time=round(timeCoefficient*i,6)

    global Coefficients, transmissionCoefficient, reflectionCoefficient, sigma, energy
    
    Function = create_function(time, axisLength, numberOfStates, Coefficients,
                                   StationaryStates, Energy, timeCoefficient)

    transmissionCoefficient, reflectionCoefficient = absorption(Function, AbsorptionFilter, absorptionLength,axisLength,
                                                                    transmissionCoefficient, reflectionCoefficient, step)


    Coefficients, coefficentsSum = create_coefficients(Function, StationaryStates, step, numberOfStates)

    functionLine.set_data(Axis,(Function*Function.conjugate()).real)


    print time, coefficentsSum + transmissionCoefficient + reflectionCoefficient


    textD = ('t = {:.2f}'.format(time) + '  T = {:.3f}'.format(transmissionCoefficient) + 
             '  R = {:.3f}'.format(reflectionCoefficient) + '  '+r'$\sigma$' +' = {:.2f}'.format(sigma) +
             '  ' +r'$\cal E$' + ' = {:.2f}'.format(energy) + r' $a$' + ' = {:.2f}'.format(potentialWidth))

    timeText.set_text(textD) 

    
    return functionLine, timeText,


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=ANIMATION_DURATION, interval=intervalTime, blit=True)

if potentialDisplay == True:
    i = np.arange(axisLength)
    plt.fill_between(Axis[i], Potential[i]*yRange, 0, alpha=0.3, color='green')
    
plt.grid(gridDisplay)

if os.path.isfile(directory + animationName + '.mp4'):
    for i in range(1,1000):
        if not os.path.isfile(directory + animationName + '(' + str(i) + ').mp4'):
            animationName = animationName + '(' + str(i) + ').mp4'
            break

else:
    animationName = animationName + '.mp4'
    
print('saving')
Writer = animation.writers['ffmpeg']
writer = Writer(fps=fpsNumber, metadata=dict(artist='Lukasz Kostrzewa'), bitrate=bitRate)
anim.save((directory + animationName), writer=writer)

print('animation saved')
exit()