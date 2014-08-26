import dictionary
import numpy as np
from brian import *

epochs = 100 
spikeInterval = 100 * ms
dictionaryLongitude = 4
spikesPerChar=3
totalTime = epochs * spikesPerChar * dictionaryLongitude * spikeInterval 

spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
LIK = SpikeGeneratorGroup(15, spiketimes)


taum = 20 * ms
taue = 1 * ms
taui = 10 * ms
Vt = 10 * mV
Vr = 0 * mV

eqs = Equations('''
      dV/dt  = (-V+ge-gi)/taum : volt
      dge/dt = -ge/taue        : volt
      dgi/dt = -gi/taui        : volt
      ''')
ADDS = NeuronGroup(N=4, model=eqs,threshold=Vt, reset=Vr)


exhitatory = Connection(LIK, ADDS , 'ge')
Wexhitatory = np.random.uniform(0,100,[15,4]) * mV
exhitatory.connect(LIK,ADDS,Wexhitatory)

inhibitory = Connection(ADDS, ADDS , 'gi')

#Connect adds layer via lateral inhibitory connections
#the diagonal should be 0 to 
Winhibitory = np.random.uniform(0,10,[4,4]) * mV
diagonal = np.diag_indices(Winhibitory.shape[0])
Winhibitory[diagonal] = 0;

inhibitory.connect(ADDS,ADDS,Winhibitory)

M = SpikeMonitor(ADDS)
Mv = StateMonitor(ADDS, 'V', record=True)
Mge = StateMonitor(ADDS, 'ge', record=True)
Mgi = StateMonitor(ADDS, 'gi', record=True)

run(totalTime)

neuronToPlot = 1
subplot(211)
raster_plot(M, title='The gupta network', newfigure=False)
subplot(223)
plot(Mv.times / ms, Mv[neuronToPlot] / mV)
xlabel('Time (ms)')
ylabel('V (mV)')
subplot(224)
plot(Mge.times / ms, Mge[neuronToPlot] / mV)
plot(Mgi.times / ms, Mgi[neuronToPlot] / mV)
xlabel('Time (ms)')
ylabel('ge and gi (mV)')
legend(('ge', 'gi'), 'upper right')
show()
