from dictionary import *
from report_results_and_logging import *
import numpy as np
from sympy import * # used for differentiation
import math # used for natural log
#from random import random # for random gen
from decimal import Decimal # Used for results reporting calculations
from brian2 import *

epochs = 100 
spikeMiliseconds = 100
spikeInterval = spikeMiliseconds * ms
spikeIntervalUnformatted = spikeMiliseconds * .001
dictionaryLongitude = 4
spikesPerChar=3
totalTime = epochs * spikesPerChar * dictionaryLongitude * spikeInterval 
epochsToPrint = [0,1,2,25,50,75]
presentResults = False
logger = False
displayAllNeuronMemPotentials = False

taum = 20 * ms
taue = 1 * ms
taui = 10 * ms
Vt = 5 * mV
Vr = -0.001 * mV

Rm = 80
## added tauS constant below ##
tauS = 2
tauMax = 30
tauMin = 2
#tauM = 30 * ms # appears different than taum above due to the article specifying this as 30
tauM = .03 # removed ms for compatibility

SimulationDuration = 10000
epochMsDuration = (SimulationDuration / epochs) * 10 # Times 10 is to adjust to Ms scale

numberOfPixels = 15
#neuronFiringThreshold = 10 * mV
neuronFiringThreshold = 10 # removed mV for compatibility
W = np.random.uniform(0.5,1.0,[4,15]) # Initial weights
# none to 1 below
numberOfNeurons = 4
R = [[1] * numberOfPixels]*dictionaryLongitude # Initial Resistance Values
testR = [[1] * numberOfPixels]*dictionaryLongitude # Initial Resistance Values
## Added 1 init below but I don't really know what the right init value is
Id = [[-.001] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
testId = [[-.001] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
## Added 1 init below but I don't really know what the right init value is
tauD = [[1] * numberOfPixels]*dictionaryLongitude # Initial tauD
testTauD = [[1] * numberOfPixels]*dictionaryLongitude # Initial tauD
tF = [[None] * numberOfPixels]*dictionaryLongitude # Initial pre-synaptic spike times

Is = [-.001] * dictionaryLongitude
Um = [-.001] * dictionaryLongitude
testUm = [-.001] * dictionaryLongitude

Ureset = -.001

APlus = 0.1
AMinus = -0.105
TauPlus = 1
TauMinus = 1
LearningRate = 0.1
ActionPotentialThreshold = .01
