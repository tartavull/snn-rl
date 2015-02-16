"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
# reference: http://matplotlib.sourceforge.net/api/animation_api.html
from brian2 import *
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import pyplot as plt
from matplotlib import animation
import h5py

class Anim3dScatterPlot:
    plottingTimeRange = 300# ms to show at a time
    def openHdf5(self, hdf5name , permision = "w"):
        self.hdf = h5py.File(hdf5name, permision)

    def closeHdf5(self):
        self.hdf.close()

    def getData(self, group,  level ):

        return self.hdf['/'+group+'/'+str(level)]

    def getDownsampleData(self, group, start, end , screenSize = 1024):
        bestLevel = 1
        data = self.getData(group,bestLevel)
        print group + " has been downsampled "+ str(bestLevel) + " times"
        start = int(start/bestLevel)
        end = int(end/bestLevel)
        if (end > data.shape[0]):
            print data.shape
            print "WARNING: you requested end = "+str(float(end))+" seconds but data length is "+str(float(data.shape[0] *bestLevel/10000))+" seconds"
            end=data.shape[0] -1

        return data[start:end]

    def incorperateData(self, group, start, end):

        #because timestep is 0.1 mseconds, and start and end is expressed in seconds.
        start = int(start * 10000)
        end = int(end * 10000)

        voltage = self.getDownsampleData('voltage', start, end)
        spikes = self.getDownsampleData('spikes', start, end)            

        return voltage, spikes

    def create3dAnim(self):
        plottingTimeRange = self.plottingTimeRange
        self.openHdf5('simulation.hdf5',permision = "r")
        start = .001 #seconds
        end = .3#.13
        voltage, spikes = self.incorperateData('voltage',start,end)
        times = np.array(range(0,len(voltage)))*.001

        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)

        # Setting the axes properties
        ax.set_xlim3d([0.0, 3.0])
        ax.set_xlabel('Neuron')
        ax.set_ylim3d([-0.03, 0.3])
        ax.set_ylabel('Time in ms*10')
        ax.set_zlim3d([-70.0, 15.0])
        ax.set_zlabel('Voltage')
        ax.set_title('Neuron voltage')
        ax.view_init(elev=15.0, azim=-10)

        c = np.array([[0],[0],[0],[0]])
        line = ax.plot(c[0, 0:1],c[1, 0:1],c[2, 0:1], lw=4)[0]
        line2 = ax.plot(c[0, 0:1],c[1, 0:1],c[2, 0:1], lw=4)[0]
        line3 = ax.plot(c[0, 0:1],c[1, 0:1],c[2, 0:1], lw=4)[0]
        line4 = ax.plot(c[0, 0:1],c[1, 0:1],c[2, 0:1], lw=4)[0]
        xVals = np.array([[0.0]*plottingTimeRange,[1]*plottingTimeRange,[2]*plottingTimeRange,[3]*plottingTimeRange])

        def generateLine(timePoint):
            timePointEnd = timePoint+plottingTimeRange
            line.set_data(xVals[0],times[0:plottingTimeRange])
            line.set_3d_properties(transpose(voltage)[0][timePoint:timePointEnd])
            line2.set_data(xVals[1],times[0:plottingTimeRange])
            line2.set_3d_properties(transpose(voltage)[1][timePoint:timePointEnd])
            line3.set_data(xVals[2],times[0:plottingTimeRange])
            line3.set_3d_properties(transpose(voltage)[2][timePoint:timePointEnd])
            line4.set_data(xVals[3],times[0:plottingTimeRange])
            line4.set_3d_properties(transpose(voltage)[3][timePoint:timePointEnd])

            if timePoint == 100 or timePoint == 200 or timePoint == 400 or timePoint == 600 or timePoint == 800:
                print 'rendering timePoint:\t',timePoint

        line_ani = animation.FuncAnimation(fig, generateLine, frames=1000, interval=.01, blit=False)
        ext = 'mp4'
        fps = 30
        codec = {'mp4': 'libx264', 'webm': 'libvpx'}.get(ext, 'mpeg4')
        line_ani.save('testVoltage.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        plt.show()

    def __init__(self):
        self.create3dAnim()
        self.closeHdf5

def main():
    run3dAnim = Anim3dScatterPlot()

if  __name__ =='__main__':main()

print 'done'