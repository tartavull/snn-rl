from dictionary import *
from report_results_and_logging import *
import numpy as np
from brian import *

epochs = 100 
spikeMiliseconds = 100
spikeInterval = spikeMiliseconds * ms
spikeIntervalUnformatted = spikeMiliseconds * .001
dictionaryLongitude = 4
spikesPerChar=3
totalTime = epochs * spikesPerChar * dictionaryLongitude * spikeInterval 
epochsToPrint = [0,1,2,25,50,75]
presentResults = True
logger = True

taum = 20 * ms
taue = 1 * ms
taui = 10 * ms
Vt = 5 * mV
Vr = 0 * mV

SimulationDuration = 10000
epochMsDuration = (SimulationDuration / epochs) * 10 # Times 10 is to adjust to Ms scale
