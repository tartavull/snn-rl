from dictionary import *
import numpy as np
from sympy import * # used for differentiation
import math # used for natural log
from decimal import Decimal # Used for results reporting calculations
from brian2 import *

epochs = 100 
spikeMiliseconds = 100
dictionaryLongitude = 26
testingSpikesPerChar=3
spikeInterval = spikeMiliseconds * ms
totalTestingTime = testingSpikesPerChar * dictionaryLongitude
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
tauM = .03 # removed ms for compatibility
SimulationDuration = 10000
epochMsDuration = (SimulationDuration / epochs) * 10 # Times 10 is to adjust to Ms scale
numberOfPixels = 15
neuronFiringThreshold = 10 # removed mV for compatibility
W = np.random.uniform(0.5,1.0,[dictionaryLongitude,15]) # Initial weights
# none to 1 below
numberOfNeurons = dictionaryLongitude
## Added 0 init below but I don't really know what the right init value is
Id = [[0] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
## Added 1 init below but I don't really know what the right init value is
tF = [[None] * numberOfPixels]*dictionaryLongitude # Initial pre-synaptic spike times
Is = [-.001] * dictionaryLongitude
Um = [-.001] * dictionaryLongitude
Ureset = -.001
APlus = 0.1
AMinus = -0.105
TauPlus = 1
TauMinus = 1
LearningRate = 0.1
ActionPotentialThreshold = .01
neuronIndex = 0
SummedDendriteGroup = 0;
SynapseToSoma = 0;
SigmaWDSyn = .0375 * mV
SigmaWDSynMod = .5
RmEq = Rm #Rm * mV
DendR = .0375 * mV
DendW = 1
DendDirac = 1	
currentEpoch = 0
totalSpikeIntervals = testingSpikesPerChar * dictionaryLongitude
IdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
IsSpikeFired = np.array([False]*dictionaryLongitude); 
IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
IsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
UmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
refractoryGraphingSwitch = False
M = None
UmM = None
weightMonitors = [None]*dictionaryLongitude
spikeIntervalCounter = 0.0
dendCalcScaling = 1.0
somaDirectScaling = 1.0
diracScaling = 1.0#1.25#0.9#1.0#0.1#0.5#1.0#0.1#1.0#0.0007#0.0007#0.0825#0.0325#0.0825
weightScaling = 1.0#1.25
lateralInhibitScaling = 0.2#0.28#0.4#0.28
lateralInhibitScaling2 = 0.7#0.8#1.0#1.5#1.7#2.0#1.8#1.3#0.6#1.6#1.1#1.6#1.3#0.7#2.5#0.5#0.8#1.5
generalClockDt = 1.0*ms#0.1*ms
positiveWeightReinforcement = 1.0#2.0#2.0#7.0#1.0#7.0#20.0#5.0#10.0#1.0#2.0
negativeWeightReinforcement = 2.0#5.0#3.5#7.5#30.0#30#2.0#30#7.5#15.0#30.0#8.0#6.0#1.0#50.0#10.0#50.0#10.0#3.0#0.5
LearningRate = 0.2#0.1#1.0#0.7
evaluateClassifier = True
accelerateTraining = False
runTimeScaling = 10
runTime = 810*ms#80*ms#810*ms#20*ms#3000*ms
#runTime = 10*ms+(10*ms*totalTestingTime) # For testing runs

class timeAndRefrac(): 	
	timeStepInterval = 0.001	
	time = 0.0000 - timeStepInterval
	refractoryPointCounter = 0.0
	refractoryPointSwitch = False
	neuronIndexCounter = 0.0
	neuronIndexSwitch = False
	lastSpikeInterval = -0.001#0.0
	spikeMiliseconds = 100
	spikeIntervalUnformatted = spikeMiliseconds * .001
	beginRefrac = [0*mV]*dictionaryLongitude	
	absRefracTime = [0]*dictionaryLongitude	
	relRefracTime = [0]*dictionaryLongitude
	relRefracTimeActivated = [False]*dictionaryLongitude	
	absRefracTimeActivated = [False]*dictionaryLongitude	
	relRefracTimeDuration = 0.03#.1#.01		
	absRefracTimeDuration = 0.02#.02#.002	
	def __init__(self):
		initalize = True	

class testRun():
	truePositiveSpikeResults = 0
	falsePositiveSpikeResults = 0
	trueNegativeSpikeResults = 0
	falseNegativeSpikeResults = 0
	firedSpikes = np.array([np.array([None])]*(totalTestingTime+3));
	'''testingSpikeWindows = np.array([None,None])
	for w in range(dictionaryLongitude+1):
		testingSpikeWindows = np.append(testingSpikeWindows, (np.array([2,3,4])+(w*3)))'''
	assignedSpikeWindows = np.array([None]*dictionaryLongitude)
	def __init__(self):
		initalize = True		

class latInhibSettings():
	voltageForInhibition = 0.0*mV
	inhibitionReduction = 1.0
	def __init__(self):
		initalize = True			