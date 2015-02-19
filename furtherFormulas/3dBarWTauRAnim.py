"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
# reference: http://matplotlib.sourceforge.net/api/animation_api.html
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pylab as pylab

tauMax = 30
tauMin = 2
neuronFiringThreshold = 10
Rm = 80
tauM = .03

class Anim3dWTauRPlot:
    def create3dBarChart(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(131, projection='3d')
        ax2 = fig.add_subplot(132, projection='3d')
        ax3 = fig.add_subplot(133, projection='3d')

        ax1.view_init(elev=10.0, azim=85)
        ax2.view_init(elev=10.0, azim=85)    
        ax3.view_init(elev=10.0, azim=85)   

        xpos = np.array([0.0])
        ypos = np.array([0.0])
        zpos = np.array([0.0])
        dx = np.array([1.0])
        dy = np.array([1.0])

        def generateLine(timePoint):
            ax1.clear()
            ax1.grid(b=False)
            ax2.clear()
            ax2.grid(b=False)
            ax3.clear()
            ax3.grid(b=False)
            ax1.set_title('Weight')
            ax2.set_title('Tau')
            ax3.set_title('Resistance')
            ax1.set_zlim([0.0, 1.0])            
            ax2.set_zlim([0.0, 30.0])            
            ax3.set_zlim([0.0, 3.5])

            w = 1-(timePoint * .01)
            tau = (tauMax - (abs(w) * (tauMax - tauMin)))
            r =  (((tau * neuronFiringThreshold) / Rm) * ((tauM / tau) ** (tauM / (tauM - tau))))

            ax1.bar3d(xpos, ypos, zpos, dx, dy, w, color=['r'], alpha=0.5)
            ax2.bar3d(xpos, ypos, zpos, dx, dy, tau, color=['g'], alpha=0.5)
            ax3.bar3d(xpos, ypos, zpos, dx, dy, r, color=['b'], alpha=0.5)

            if timePoint == 10 or timePoint == 20 or timePoint == 40 or timePoint == 60 or timePoint == 80:
                print 'rendering timePoint:\t',timePoint

        line_ani = animation.FuncAnimation(fig, generateLine, frames=100, interval=1, blit=False)        
        ext = 'mp4'
        fps = 30
        codec = {'mp4': 'libx264', 'webm': 'libvpx'}.get(ext, 'mpeg4')
        line_ani.save('WeightTauRes.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

        plt.show()

    def __init__(self):
        self.create3dBarChart()

def main():
    run3dAnim = Anim3dWTauRPlot()

if  __name__ =='__main__':main()

print 'done'