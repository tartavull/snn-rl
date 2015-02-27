#Just to know the last time this was executed
import time
print time.ctime()
import h5py
from bisect import bisect_left
import matplotlib.pyplot as plt
from brianPlotter import *
from gupta_paper_further_formulas_brianversion19 import *

# ***                                               ***
# Be sure to set correct brian main program import file
# ***                                               ***

bp = brianPlotter('simulation.hdf5')

gupta_paper = gupta_paper()

# spikeVoltage and testSpikeVoltage used for reformatting to avoid formatting issue with plotting
spikeVoltage = []
for timeIndex in range(len(gupta_paper.M.i[:])):
	spikeVoltage.append(gupta_paper.M.i[timeIndex].tolist())

spikes = zip(spikeVoltage, gupta_paper.M.t)
#testSpikes = zip(testSpikeVoltage, gupta_paper.testM.t)

bp.saveData('spikes', np.asarray(spikes), downsample=False);
#bp.saveData('testSpikes', np.asarray(testSpikes), downsample=False);

voltage =  np.asarray(gupta_paper.UmM.v2.T/mV)
#testVoltage =  np.asarray(gupta_paper.testUmM3.v2.T/mV)
#print 'voltage\t',voltage
bp.saveData('voltage',voltage)
#bp.saveData('testVoltage',testVoltage)
for weightIndex in range(dictionaryLongitude):
	bp.saveData('weights'+str(weightIndex),np.asarray(gupta_paper.weightMonitors[weightIndex].w.T/volt),downsample=False)
'''weights0 =  np.asarray(gupta_paper.WM0.w.T/volt)
weights1 =  np.asarray(gupta_paper.WM0.w.T/volt)
weights2 =  np.asarray(gupta_paper.WM0.w.T/volt)
weights3 =  np.asarray(gupta_paper.WM0.w.T/volt)
bp.saveData('weights0',weights0,downsample=False)
bp.saveData('weights1',weights1,downsample=False)
bp.saveData('weights2',weights2,downsample=False)
bp.saveData('weights3',weights3,downsample=False)'''

del bp #you have to delete the object so the file is closed

print time.ctime()