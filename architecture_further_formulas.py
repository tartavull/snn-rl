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
#spikeIntervalUnformatted = spikeMiliseconds * .001
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
#Id = [[-.001] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
Id = [[0] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
testId = [[-.001] * numberOfPixels]*dictionaryLongitude # Initial Dendritic Post Synaptic Current
## Added 1 init below but I don't really know what the right init value is
tauD = [[1] * numberOfPixels]*dictionaryLongitude # Initial tauD
testTauD = [[1] * numberOfPixels]*dictionaryLongitude # Initial tauD
tF = [[None] * numberOfPixels]*dictionaryLongitude # Initial pre-synaptic spike times

Is = [-.001] * dictionaryLongitude
Um = [-.001] * dictionaryLongitude
testUm = [-.001] * dictionaryLongitude

Ureset = -.001
#Ureset = 0

APlus = 0.1
AMinus = -0.105
TauPlus = 1
TauMinus = 1
LearningRate = 0.1
ActionPotentialThreshold = .01

#absRefracTime = [0]*dictionaryLongitude	
#relRefracTime = [0]*dictionaryLongitude
#timeStepInterval = 0.001	
#time = 0.0000 - timeStepInterval
neuronIndex = 0
#lastSpikeInterval = -0.001#0.0
SummedDendriteGroup = 0;
SynapseToSoma = 0;
SigmaWDSyn = .0375 * mV
SigmaWDSynMod = .5
RmEq = Rm #Rm * mV
DendR = .0375 * mV
DendW = 1
DendDirac = 1	
currentEpoch = 0
'''refractoryPointCounter = 0.0
testRefractoryPointCounter = 0.0
refractoryPointSwitch = False
testRefractoryPointSwitch = False
neuronIndexCounter = 0.0
neuronIndexSwitch = False'''
totalSpikeIntervals = 12	
testSpikesFiredInInterval = np.array([[False]*(totalSpikeIntervals+1)]*dictionaryLongitude)
IdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
testIdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
IsSpikeFired = np.array([False]*dictionaryLongitude); 
testIsSpikeFired = np.array([False]*dictionaryLongitude);	
IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
testIdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
IsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
testIsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
UmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
testUmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
refractoryGraphingSwitch = False

M = None
testM = None
UmM3 = None
testUmM3 = None
spikeIntervalCounter = 0.0
#relRefracTimeActivated = [False]*dictionaryLongitude	
#absRefracTimeActivated = [False]*dictionaryLongitude
# Parameters for tuning
#relRefracTimeDuration = 0.03#.1#.01		
#absRefracTimeDuration = 0.02#.02#.002	
#dendCalcScaling = 0.16#.32#1.0
dendCalcScaling = 1.0
#somaDirectScaling = .0005#.0005#.001#.0005#.001#0.002285#0.002857
somaDirectScaling = 1.0
diracScaling = 1.0#1.25#0.9#1.0#0.1#0.5#1.0#0.1#1.0#0.0007#0.0007#0.0825#0.0325#0.0825
weightScaling = 1.0#1.25
#diracScaling = 1.0
#uMScaling = 1.0#.32
lateralInhibitScaling = 0.2#0.28#0.4#0.28
lateralInhibitScaling2 = 0.7#0.8#1.0#1.5#1.7#2.0#1.8#1.3#0.6#1.6#1.1#1.6#1.3#0.7#2.5#0.5#0.8#1.5
generalClockDt = 1.0*ms#0.1*ms
positiveWeightReinforcement = 1.0#2.0#2.0#7.0#1.0#7.0#20.0#5.0#10.0#1.0#2.0
negativeWeightReinforcement = 2.0#5.0#3.5#7.5#30.0#30#2.0#30#7.5#15.0#30.0#8.0#6.0#1.0#50.0#10.0#50.0#10.0#3.0#0.5
LearningRate = 0.2#0.1#1.0#0.7
evaluateClassifier = True
accelerateTraining = False
runTime = 132*ms#1000*ms#300*ms#400*ms#155*ms#42*ms#135*ms#72*ms#135*ms#42*ms#100*ms#500*ms#5000*ms#100*ms#22*ms#15*ms#135*ms#15*ms#130*ms

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
	'''timeStepInterval = self.timeStepInterval
	time = self.time
	refractoryPointCounter = self.refractoryPointCounter
	refractoryPointSwitch = self.refractoryPointSwitch
	neuronIndexCounter = self.neuronIndexCounter
	lastSpikeInterval = self.lastSpikeInterval
	spikeIntervalUnformatted = self.spikeIntervalUnformatted'''
	def __init__(self):
		initalize = True	

#timeAndRefrac = timeAndRefrac()

class testRun():
	truePositiveSpikeResults = 0
	falsePositiveSpikeResults = 0
	trueNegativeSpikeResults = 0
	falseNegativeSpikeResults = 0
	correctSpikes = np.array([[1]*(totalSpikeIntervals+1)]*dictionaryLongitude);#[[1] * totalSpikeIntervals]*dictionaryLongitude
	correctOrderOfNeurons = [2,0,1,3]
	correctSpikes[correctOrderOfNeurons[0]][0] = 0
	correctSpikes[correctOrderOfNeurons[0]][4:13] = 0
	correctSpikes[correctOrderOfNeurons[1]][0:4] = 0
	correctSpikes[correctOrderOfNeurons[1]][7:13] = 0
	correctSpikes[correctOrderOfNeurons[2]][0:7] = 0
	correctSpikes[correctOrderOfNeurons[2]][10:13] = 0
	correctSpikes[correctOrderOfNeurons[3]][0:10] = 0
	def __init__(self):
		initalize = True		

class latInhibSettings():
	voltageForInhibition = 0.0*mV
	inhibitionReduction = 1.0
	def __init__(self):
		initalize = True			
