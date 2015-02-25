import h5py
import numpy as np

class processHDF5Data:
    def __init__(self, hdf5name , permision = "w"):
        self.hdf = h5py.File(hdf5name, permision)

    def __del__ (self):
        self.hdf.close()

    def getData(self, group,  level ):
        return self.hdf['/'+group+'/'+str(level)]

    def getDownsampleData(self, group, start, end, screenSize = 1024):
        bestLevel = 1
        data = self.getData(group,bestLevel)
        print group + " has been downsampled "+ str(bestLevel) + " times"
        start = int(start/bestLevel)
        end = int(end/bestLevel)
        if (end > data.shape[0]):
            print data.shape
            print "WARNING: you requested end = "+str(float(end))+" seconds but data length is "+str(float(data.shape[0] *bestLevel/10000))+" seconds"
            end=data.shape[0] -1
            print 'end',end

        return data[start:end]

    def incorperateData(self, timeSeriesData, names, start, end):
        #because timestep is 0.1 mseconds, and start and end is expressed in seconds.
        start = int(start * 10000)
        end = int(end * 10000)

        for inputDataIndex in range(len(timeSeriesData)):
            timeSeriesData[inputDataIndex] = self.getDownsampleData(names[inputDataIndex], start, end)

        return timeSeriesData