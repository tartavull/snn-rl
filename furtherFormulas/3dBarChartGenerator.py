# from http://pythonprogramming.net/3d-bar-charts-python-matplotlib/
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pylab as pylab
 
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
 
xpos = np.ones(60)
xpos[15:30] *= 2
xpos[30:45] *= 3
xpos[45:60] *= 4
ypos = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

num_elements = len(xpos)
zpos = np.zeros(60)
dx = np.ones(60)*.1
dy = np.ones(60)*.5
dz = [0.90433194, 0.6139531, 0.50387484, 0.55220372, 0.51213536, 0.85443374, 0.99955922, 0.5039825, 0.73091913, 0.9780236, 0.5241028, 0.71571812, 0.93782861, 0.51210244, 0.73074697,
	0., 0., 0.03412608, 0., 0.90455366, 0.78683668, 0., 0.95912629, 0.7282637, 0., 0.78548583, 0.78935491, 0.03193823, 0.00609877, 0.17287094,
	0.4474444, 0., 0.98135641, 0., 0.96315942, 0., 0., 0., 0.15930208, 0., 0.77299245, 0., 0., 0.71739497, 0.02804206,
	0., 0., 0.99815102, 0., 0.9239562, 0., 0., 0.32862838, 0.29682383, 0., 0.85108903, 0., 0., 0., 0.6687179]

colors = ['r']*15+['g']*15+['b']*15+['y']*15

xLabel = ax1.set_xlabel('\nOutput Neuron', linespacing=1.2)
yLabel = ax1.set_ylabel('\nInput Neuron', linespacing=1.2)
zLabel = ax1.set_zlabel('\nWeight', linespacing=1.2)

neuron1 = plt.Rectangle((0, 0), 0.1, 0.1,fc='r')
neuron2 = plt.Rectangle((0, 0), 0.1, 0.1,fc='g')
neuron3 = plt.Rectangle((0, 0), 0.1, 0.1,fc='b')
neuron4 = plt.Rectangle((0, 0), 0.1, 0.1,fc='y')

ax1.legend((neuron1,neuron2,neuron3,neuron4),("neuron 1","neuron 2","neuron 3","neuron 4"),'best')

ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, alpha=0.5)
ax1.view_init(elev=60.0, azim=40)

pylab.savefig('Weights3dBar.jpg')

plt.show()
