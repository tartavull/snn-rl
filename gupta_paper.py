import dictionary
import numpy as np
from brian import *

spiketimes = dictionary.spikeTimes(dictionaryLongitude=4, spikeInterval=100 * ms, spikesPerChar=3, epochs=100)
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
Winhibitory = np.random.uniform(0,10,[4,4]) * mV
inhibitory.connect(ADDS,ADDS,Winhibitory)


Mv = StateMonitor(ADDS, 'V', record=True)
Mge = StateMonitor(ADDS, 'ge', record=True)
Mgi = StateMonitor(ADDS, 'gi', record=True)

run(1000 * ms)

figure()
subplot(211)
plot(Mv.times / ms, Mv[0] / mV)
subplot(212)
plot(Mge.times / ms, Mge[0] / mV)
plot(Mgi.times / ms, Mgi[0] / mV)
show()
