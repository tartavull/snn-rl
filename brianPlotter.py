import h5py
from bisect import bisect_left
import matplotlib.pyplot as plt
import numpy as np

class brianPlotter:

    def __init__(self, hdf5name , permision = "w"):
        self.hdf = h5py.File(hdf5name, permision)

    def saveData(self, dataName , dataArray , downsample=True):
        self.grp = self.hdf.create_group(dataName)
        self.grp.create_dataset("1", data=dataArray)
        
        if (downsample):
            level = 1
            step = 10
            while ( dataArray.shape[0]/level > 1000 ):

                self.minMaxDownsampling(dataName, level , step)
                level = level * step

        print dataName + ' sucessfully saved';
    def __del__ (self):
        self.hdf.close()
        
    def minMaxDownsampling(self, group, baseLevel , step):
        data = self.hdf['/'+group+'/'+str(baseLevel)]
        
        shape = list(data.shape)
        shape[0] = shape[0] / step
        downsample = np.empty(shape) 
        
        for bin in range(0,data.shape[0],step):
            index = bin/step
            if (index % 2 == 0):
                downsample[index,:] = data[bin:(bin+step)].max(axis=0)
            else:  
                downsample[index,:] = data[bin:(bin+step)].min(axis=0)
        level = int(baseLevel)*step
        self.grp.create_dataset(str(level), data=downsample)
    
    def plotLine(self, group, start, end):

        data = self.getDownsampleData(group, start, end)
        length = self.getGroupLength(group)
        
        nrows = data.shape[1]
     
        fig, axes = plt.subplots(ncols=1, nrows=nrows, figsize=(16,nrows*2.5),sharey=True,sharex = True)
        fig.suptitle(group, fontsize=20)

        for index in range(nrows):     
            legend = 'neuron '+str(index)
            if (len(data.shape) == 3):
                legend = []
                for l in range(data.shape[2]):
                    legend.append('Synapse '+str(index)+ '-'+str(l))
            
            axes[index].plot(data[:,index])
            axes[index].legend(legend,framealpha=0.5)

    def plotScatter(self, start, end):
        
        spikes = self.getData('spikes',1)

        start =  self.binary_search(spikes[:,1], start)
        end = self.binary_search(spikes[:,1], end)

        fig, ax = plt.subplots(subplot_kw=dict(axisbg='#EEEEEE'),figsize=(16,4))
        points = plt.scatter(spikes[start:end,1],spikes[start:end,0], s=100)
        plt.yticks(np.arange(0,4, 1.0))
        plt.ylabel("Neuron number")
        plt.xlabel("Time(seconds)")
        plt.grid(True)

    def getData(self, group,  level ):

        return self.hdf['/'+group+'/'+str(level)]

    def getDownsampleData(self, group, start, end , screenSize = 1024):
       
        bestLevel = self.getBestLevel(group, start , end, screenSize)
        
        data = self.getData(group,bestLevel)
        print "data has been downsampled "+ str(bestLevel) + " times"
        return data[start:end]
            
    def getLevels(self, group):

        group = self.hdf['/'+group]
        ilist = group.items()
        levels = list()

        for index in range(len(ilist)):
            levels.append(int(ilist[index][0]))
        
        return levels
    
    def getBestLevel(self, group, start, end , screenSize = 1024):
        
        requireLevel = (end - start) / screenSize;
        possibleLevels = self.getLevels(group)
        
        min = abs(requireLevel - possibleLevels[0])
        bestFit = 0
        for i in range(1,len(possibleLevels)):
            if (min > abs(requireLevel - possibleLevels[i])):
                min = abs(requireLevel - possibleLevels[i])
                bestFit = i        
       
        
        return possibleLevels[bestFit]
    
    def getGroupLength(self,group): #not yet implemented
        return 1200000
    
    def binary_search(self, a, x, lo=0, hi=None): 
        x = x/1000
        hi = hi if hi is not None else len(a)  
        return bisect_left(a,x,lo,hi) 