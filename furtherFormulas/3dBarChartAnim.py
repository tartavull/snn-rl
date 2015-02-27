"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
# reference: http://matplotlib.sourceforge.net/api/animation_api.html
# http://pythonprogramming.net/3d-bar-charts-python-matplotlib/
from brian2 import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab as pylab
from decimal import Decimal
from gupta_paper_further_formulas_brianversion19 import *
from processHDF5Data import *

plottingTimeRange = 300# ms to show at a time

class Anim3dBarChart:
    plottingTimeRange = self.plottingTimeRange

    def create3dBarChart(self):
        plottingTimeRange = self.plottingTimeRange
        hdf5Data = processHDF5Data('../simulation.hdf5',permision = "r")
        start = .001 #seconds
        end = 5.0#.5#.13
        weightMonitors = [None]*dictionaryLongitude
        names = [None]*dictionaryLongitude
        for weightIndex in dictionaryLongitude:
            names[weightIndex] = 'weights'+str(weightIndex)
        weightMonitor = hdf5Data.incorperateData(weightMonitors,names,start,end)
        times = np.array(range(0,len(weightMonitor)))*.001

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')
         
        xpos = np.ones(60)
        xpos[0:15] *= 0
        xpos[15:30] *= 1
        xpos[30:45] *= 2
        xpos[45:60] *= 3
        ypos = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

        num_elements = len(xpos)
        zpos = np.zeros(60)
        dx = np.ones(60)*.1
        dy = np.ones(60)*.5

        colors = ['r']*15+['g']*15+['b']*15+['y']*15

        ax1.set_xlim([0.0, 3.0])
        ax1.set_xlabel('Neuron')
        ax1.set_ylim([0, 16.0])
        ax1.set_ylabel('Time in ms*10')
        ax1.set_zlim([0.0, 1.0])
        ax1.set_zlabel('Voltage')
        ax1.set_title('Neuron voltage')
        ax1.view_init(elev=60.0, azim=40)

        self.weightTotalTime = 0

        def generateLine(timePoint):
            self.weightTotalTime += 50
            ax1.grid(b=False)
            ax1.clear()
            ax1.grid(b=False)
            ax1.set_xlim([0.0, 3.0])
            ax1.set_ylim([0, 16.0])
            ax1.set_zlim([0.0, 1.0])
            ax1.spines['bottom'].set_color('#dddddd')
            ax1.spines['top'].set_color('#dddddd') 
            ax1.spines['right'].set_color('red')
            ax1.spines['left'].set_color('red')
            weightsCombined = []
            weightsCombined.extend(weights0[self.weightTotalTime])
            weightsCombined.extend(weights1[self.weightTotalTime])
            weightsCombined.extend(weights2[self.weightTotalTime])
            weightsCombined.extend(weights3[self.weightTotalTime])            
            ax1.bar3d(xpos, ypos, zpos, dx, dy, weightsCombined, color=colors, alpha=0.5)

        line_ani = animation.FuncAnimation(fig, generateLine, frames=90, interval=1, blit=False)        
        ext = 'mp4'
        fps = 30
        codec = {'mp4': 'libx264', 'webm': 'libvpx'}.get(ext, 'mpeg4')
        line_ani.save('trainingWeights.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        plt.show()

        del(hdf5Data)

    def __init__(self):
        self.create3dBarChart()

def main():
    run3dAnim = Anim3dBarChart()

if  __name__ =='__main__':main()

print 'done'