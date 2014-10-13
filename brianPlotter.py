import h5py
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
    
    def plot(self, group,  level ):

        data = self.hdf['/'+group+'/'+str(level)][0,:]
        time = range(0,data.shape[0]*level,level)
        plt.plot(time,data)

    def getData(self, group,  level ):

        return self.hdf['/'+group+'/'+str(level)]

    def getDownsampleData(self, group, start, end):

        screenSize = 1024
        requireLevel = (start - end) / screenSize;
        
        return false #not yet implemented

    def getLevels(self, group):

        group = self.hdf['/'+group]
        return group.name
