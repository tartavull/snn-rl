from dictionary import *
from report_results_and_logging import *
import numpy as np
from sympy import * # used for differentiation
import math # used for natural log
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

Rm = 80
tauMax = 30
tauMin = 2
tauM = 30 * ms # appears different than taum above due to the article specifying this as 30

SimulationDuration = 10000
epochMsDuration = (SimulationDuration / epochs) * 10 # Times 10 is to adjust to Ms scale

numberOfPixels = 15
neuronFiringThreshold = 10 * mV
W = np.random.uniform(0.5,1.0,[15,4]) # Initial weights
R = [[None] * dictionaryLongitude for i in range(numberOfPixels)] # Initial Resistance Values
Id = [[None] * dictionaryLongitude for i in range(numberOfPixels)] # Initial Dendritic Post Synaptic Current
tauD = [[None] * dictionaryLongitude for i in range(numberOfPixels)] # Initial tauD
tF = [[None] * dictionaryLongitude for i in range(numberOfPixels)] # Initial pre-synaptic spike times
