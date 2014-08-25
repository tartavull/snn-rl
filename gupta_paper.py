import dictionary

from brian import *

tau = 20 * msecond        # membrane time constant
Vt = -50 * mvolt          # spike threshold
Vr = -60 * mvolt          # reset value
El = -49 * mvolt          # resting potential (same as the reset)
psp = 0.5 * mvolt         # postsynaptic potential size


spiketimes = dictionary.spikeTimes(dictionaryLongitude=4, spikeInterval=100 * ms, spikesPerChar=3, epochs=100)

LIK = SpikeGeneratorGroup(15, spiketimes)
ADDS = NeuronGroup(N=4, model='dV/dt = -(V-El)/tau : volt',threshold=Vt, reset=Vr)

C1 = Connection(LIK, ADDS)

run(100 * ms)
