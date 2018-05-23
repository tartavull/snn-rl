# from http://pythonprogramming.net/3d-bar-charts-python-matplotlib/
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import animation
import numpy as np
import pylab as pylab
import math as math
from decimal import *
from savedWeights import *
 
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
colorNumber = -1
ypos = []
yposNeuronOrder = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
weights = retreiveTestW()
colors = []
neurons = [None]*len(weights)
numberOfWeights = len(weights)
numberOfPixels = 15
xpos = np.ones(numberOfPixels*numberOfWeights)
num_elements = len(xpos)
zpos = np.zeros(numberOfPixels*numberOfWeights)
dx = np.ones(numberOfPixels*numberOfWeights)*.9
dy = np.ones(numberOfPixels*numberOfWeights)*.5
timeDelay = 24#3

for weightIndex in range(numberOfWeights):
	colorNumber += 1
	if colorNumber == 4: colorNumber = 0
	if colorNumber == 0: colors = np.append(colors, ['r']*15); neurons[weightIndex] = plt.Rectangle((0, 0), 0.1, 0.1,fc='r');
	if colorNumber == 1: colors = np.append(colors, ['g']*15); neurons[weightIndex] = plt.Rectangle((0, 0), 0.1, 0.1,fc='g');
	if colorNumber == 2: colors = np.append(colors, ['b']*15); neurons[weightIndex] = plt.Rectangle((0, 0), 0.1, 0.1,fc='b');
	if colorNumber == 3: colors = np.append(colors, ['y']*15); neurons[weightIndex] = plt.Rectangle((0, 0), 0.1, 0.1,fc='y');
	ypos = np.append(ypos, yposNeuronOrder)
	xpos[(15*weightIndex):(15*(weightIndex+1))] *= (weightIndex+1)

colors = colors.tolist();ypos = ypos.tolist();xpos = xpos.tolist();

def generateGraph(timePoint):
	timePointScaled =  math.floor(Decimal(format((timePoint), '.1f'))/Decimal(format(timeDelay, '.1f')))
	ax1.clear()
	ax1.grid(b=False)	
	ax1.set_axisbelow(True)
	dz = []
	for weightIndex in range(numberOfWeights):
		if weightIndex == timePointScaled: dz = np.append(dz, weights[weightIndex])
		else: dz = np.append(dz, [0]*numberOfPixels)
	dz = dz.tolist()
	xLabel = ax1.set_xlabel('\nOutput Neuron', linespacing=1.2)
	yLabel = ax1.set_ylabel('\nInput Neuron', linespacing=1.2)
	zLabel = ax1.set_zlabel('\nWeight', linespacing=1.2)
	ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.95)
	if timePointScaled==5 or timePointScaled==10 or timePointScaled==15 or timePointScaled==20 or timePointScaled==25: print 'processing timepoint:\t',timePointScaled

ax1.view_init(elev=30.0, azim=40)

line_ani = animation.FuncAnimation(fig, generateGraph, frames=(numberOfWeights*timeDelay), interval=1, blit=False)        
ext = 'mp4'
fps = 30
codec = {'mp4': 'libx264', 'webm': 'libvpx'}.get(ext, 'mpeg4')
line_ani.save('26CharTrainedWeights.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
print 'Completed saving the video'

#plt.show()
