'''
	Copyright 2015, Nate Sutton and Ignacio Tartavull
	This is the main file for the Spiking Neual Networks
	Reinforcement Learning simulation.  

	More info:
	https://github.com/tartavull/snn-rl/blob/master/README.md
	http://nbviewer.ipython.org/github/tartavull/snn-rl/blob/master/notebooks/introduction.ipynb
	http://nbviewer.ipython.org/github/tartavull/snn-rl/blob/master/FFSSN.ipynb
'''

from furtherFormulas.architecture_further_formulas import *
from furtherFormulas.cofactorCalculations import *
from furtherFormulas.timeAndRefracCalcs import *
from furtherFormulas.outputPrinting import *
from furtherFormulas.testingProcesses import *
from furtherFormulas.lateralInhibition import *
from furtherFormulas.generalUtilities import *

timeAndRefrac = timeAndRefrac	

class gupta_paper:
	'''
		Main program variables are set the simulation is run.  Specific neuron models are
		defined for computing processing of electrophysiology in the soma, direct to soma
		signals, and active dendrites.  Equations with variables
		are defined in the equations() objects.
	'''
	neuralnet = Network()
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, testingSpikesPerChar, testingEpochs)
	trainingSpikeTimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, trainingSpikesPerChar, trainingEpochs)
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	# W = W and other lines below are to avoid an odd 'reference before assignment' error
	W = W
	tauM = tauM
	neuronIndex = neuronIndex
	generalClockDt = generalClockDt
	runTime = runTime
	runTimeScaling = runTimeScaling
	evaluateClassifier = evaluateClassifier
	accelerateTraining = accelerateTraining
	diracScaling = diracScaling
	somaDirectScaling = somaDirectScaling
	negativeWeightReinforcement = negativeWeightReinforcement
	positiveWeightReinforcement = positiveWeightReinforcement
	timeAndRefrac = timeAndRefrac	
	testRun = testRun
	latInhibSettings = latInhibSettings
	standardPrint = standardPrint
	verbosePrint = verbosePrint
	testingRunTime = testingRunTime
	optResultsFile = optResultsFile
	minWeightsRand = minWeightsRand
	maxWeightsRand = maxWeightsRand
	totalSpikeIntervals = totalSpikeIntervals

	def run_model(self):
		neuralnet = self.neuralnet
		dictionary = self.dictionary

		eqs = Equations('''
			dv/dt = v/(1*second): volt
			dprelimV/dt = (-prelimV+((Rm/mV)*(SynI+DendI*1.0)))/(tauM) : volt (unless refractory)
			Rm = 80*mV : volt
			tauM = 30*ms : second
	        V : volt
	        DendI : volt
	        SynI : volt
	        v2 : volt	
			UmSpikeFired : volt	
			beginRefrac : volt
		    ''')			

		dendriteEqs = Equations('''
			dv/dt = (((-v/mV)+((r/mV)*(w/volt)*(dirac/volt)))/(tau))*mV : volt
			V : volt
	        r : volt
	        w : volt
	        dirac : volt
	        tau : second
	        v2: volt
			''')

		directToSomaEqs = Equations('''
			dv/dt = (((-v/mV)+(summedWandDirac/volt))/(tauS))*mV : volt
			tauS = 2*ms : second
			V : volt
			summedWandDirac : volt
			v2 : volt
			spikeFired : boolean
			''')		

		class ADDSNeuronModel(NeuronGroup, gupta_paper): 
			'''
				This is the model used for electrophysiology occuring in the Soma
			'''
			neuronIndex = self.neuronIndex
			generalClockDt = self.generalClockDt

			def __init__(self, params):
				self = parseArgs(self, sys.argv, dictionaryLongitude)		

				NeuronGroup.__init__(self, N=dictionaryLongitude, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; v2=10*mV;UmSpikeFired=1*mV;beginRefrac=1*mV;inhibitionVoltage=prelimV',refractory=8*ms,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					neuronIndex = self.neuronIndex
					timeAndRefrac.spikeIntervalCounter = (floor(timeAndRefrac.time/timeAndRefrac.spikeIntervalUnformatted) * timeAndRefrac.spikeIntervalUnformatted)*10

					def dendCalcs(neuronIndex, ADDSObj, dendObj, spiketimes, evaluationActive):
						'''
							Below sequentially Dirac, Tau, then Resistance are calculated every end of a spike-time interval.
							The resulting Dend I is added to the Um calc for the ADDS soma.
						'''
						timeAndRefrac = self.timeAndRefrac

						# Dirac
						dendObj[neuronIndex].dirac = diracCalc(dendObj, neuronIndex, spiketimes, timeAndRefrac.time, timeAndRefrac.lastSpikeInterval)

						# Initialize weights
						if (evaluationActive==False and timeAndRefrac.time == 0.000):
							dend[neuronIndex].w = W[neuronIndex]*volt

						if (evaluationActive==False and timeAndRefrac.refractoryPointSwitch==True):
							# Only change weights of neuron fired
							if ADDSObj.UmSpikeFired[neuronIndex] == 1*mV:
								# Weights
								WeightChangeCalculation(neuronIndex, spiketimes, timeAndRefrac.time, self.negativeWeightReinforcement, self.positiveWeightReinforcement, M, dendObj)	
							# Tau
							tauDCalc(neuronIndex, dendObj)
							# Resistance
							resistanceCalc(neuronIndex, dendObj, self.tauM)

						# TODO: check do I need additional loop below?
						for indexOfDend in range(dictionaryLongitude):
							ADDSObj.DendI[indexOfDend] = sum(dendObj[indexOfDend].v[:])*dendCalcScaling

						#print 'ADDSObj.t',ADDSObj.t,'ADDSObj.DendI',ADDSObj.DendI,'neuronIndex',neuronIndex,'dendObj[neuronIndex].dirac',dendObj[neuronIndex].dirac,'dendObj[neuronIndex].tau',dendObj[neuronIndex].tau,'dendObj[neuronIndex].w',dendObj[neuronIndex].w,'dendObj[neuronIndex].r',dendObj[neuronIndex].r

					def somaDirectCalcs(neuronIndex, ADDSObj, somaDirectObj, dendObj):
						dotProductWandDirac =  sum(w*d for w,d in zip(dendObj[neuronIndex].w[:], dendObj[neuronIndex].dirac[:]))
						#somaDirectObj.summedWandDirac[neuronIndex] = ((dotProductWandDirac*volt)/(volt**2))*self.somaDirectScaling
						somaDirect.summedWandDirac[neuronIndex] = ((dotProductWandDirac*volt)/(volt**2))*self.somaDirectScaling

						for neuronNumber in range(dictionaryLongitude):
							#ADDSObj.SynI[neuronNumber] = somaDirectObj.v[neuronNumber]		
							ADDS.SynI[neuronNumber] = somaDirect.v[neuronNumber]		

						#print 'ADDSObj.t',ADDSObj.t,'ADDSObj.SynI',ADDSObj.SynI,'neuronIndex',neuronIndex,'somaDirectObj.summedWandDirac',somaDirectObj.summedWandDirac[neuronIndex],'dendObj[neuronIndex].w',dendObj[neuronIndex].w,'dendObj[neuronIndex].dirac',dendObj[neuronIndex].dirac

					def mainSimulationCalcs(ADDSObj, dendObj, somaDirectObj, spiketimes, evaluationActive):
						'''
							dend then somaDirect calcs are done which are then used to set lat inhib.
							Soma Um calcs are done automatically using equations entered for brian
							once dend and somaDirect are updated
						'''
						preTNorm = self.timeAndRefrac.time
						tNorm = preTNorm - (floor((preTNorm/.001)*.01) * .1)
						
						self.timeAndRefrac = timePeriodAndRefractoryCalcs(self.timeAndRefrac)

						if (evaluationActive==True) and (timeAndRefrac.time == 0.000 or timeAndRefrac.time == 0.001):
							initializeTrainedModelParameters(dendObj)

						# Option to accelerate computations for training
						if self.accelerateTraining == False or (evaluationActive == False and (tNorm <= .005 or tNorm >= .096)):													
							if self.accelerateTraining == True and (tNorm >= .096 and tNorm < .099):
								for i in range(dictionaryLongitude):
									ADDSObj.DendI[i]=0*mV
									ADDSObj.SynI[i]=0*mV
									ADDSObj.prelimV[i]=0*mV
									ADDSObj.v[i]=0*mV								
									for i2 in range(len(dend[0].v)):
										dendObj[i].v[i2] = 0*mV
									somaDirectObj.v[i] = 0*mV		

							for neuronIndex in range(dictionaryLongitude):
								dendCalcs(neuronIndex, ADDSObj, dendObj, spiketimes, evaluationActive)

								somaDirectCalcs(neuronIndex, ADDSObj, somaDirectObj, dendObj)								

							ADDSObj.v, self.latInhibSettings = lateralInhibition(ADDSObj, self.timeAndRefrac, self.latInhibSettings)

							#if ADDSObj.t > 100*ms:
								#for i in range(dictionaryLongitude):
								#	somaDirectObj.summedWandDirac[i] += 20*mV
								#print 'ADDS.t',ADDS.t,'ADDS.SynI',ADDS.SynI,'ADDS.DendI',ADDS.DendI,'ADDSObj.v',ADDSObj.v,'somaDirectObj.summedWandDirac',somaDirectObj.summedWandDirac,'somaDirectObj.v',somaDirectObj.v

							for neuronIndex in range(dictionaryLongitude): 
								ADDSObj.v2, somaDirectObj.v2, self.timeAndRefrac = checkForResets(neuronIndex, ADDSObj, dendObj, somaDirectObj, self.timeAndRefrac)

							ADDSObj.UmSpikeFired, self.testRun = evaluateClassifierPerf2(ADDSObj, self.testRun)

						#roundedSecondsTime = math.floor(Decimal(format((ADDSObj.t), '.1f'))/Decimal(format((1.0*second), '.1f')))
						#if printWeights < roundedSecondsTime:
						#	self.printWeights = roundedSecondsTime; printWeights(dendObj);
					if self.evaluateClassifier == False:
						mainSimulationCalcs(ADDS, dend, somaDirect, self.trainingSpikeTimes, False)
					else:
						mainSimulationCalcs(ADDS, dend, somaDirect, self.spiketimes, True)					

				self.contained_objects.append(additionToNetwork)				

		class DendriteNeuronModel(NeuronGroup):
			generalClockDt = self.generalClockDt
			def __init__(self, params): 
				NeuronGroup.__init__(self, N=15, model=dendriteEqs,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					placeHolderForLaterContent = True
				self.contained_objects.append(additionToNetwork)

		class SomaDirectNeuronModel(NeuronGroup): 
			generalClockDt = self.generalClockDt
			def __init__(self, params): 
				NeuronGroup.__init__(self, N=dictionaryLongitude, model=directToSomaEqs,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					placeHolderForLaterContent = True
				self.contained_objects.append(additionToNetwork)		

		dend = [None]*dictionaryLongitude
		weightMonitors = [None]*dictionaryLongitude
		for firstLayerIndex in range(dictionaryLongitude):
			dend[firstLayerIndex] = DendriteNeuronModel(15)	
			weightMonitors[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'w', record=True)
			neuralnet.add(dend[firstLayerIndex])
			neuralnet.add(weightMonitors[firstLayerIndex])
		somaDirect = SomaDirectNeuronModel(dictionaryLongitude)
		neuralnet.add(somaDirect)
		ADDS = ADDSNeuronModel(self)			
		M = SpikeMonitor(ADDS)
		Mv = StateMonitor(ADDS, 'V', record=True)
		UmM = StateMonitor(ADDS, 'v2', record=True)
		self.M = M # for ipython compatibility
		self.UmM = UmM 
		self.weightMonitors = weightMonitors

		neuralnet.add(ADDS)
		neuralnet.add(M)
		neuralnet.add(UmM)
		if (ADDS.evaluateClassifier==True): ADDS.runTime = ADDS.testingRunTime
		ADDS.runTime *= ADDS.runTimeScaling # scaling factor
		if ADDS.standardPrint: neuralnet.run(ADDS.runTime,report='text')
		else: neuralnet.run(ADDS.runTime,report='stderr')

		OutputEvaluationResults(dend, self.testRun, ADDS.verbosePrint, ADDS.evaluateClassifier)
		accuracyPerc = totalCorrectPercentage()
		precisionPerc = precisionPercentage()
		writeOptimizationResults(ADDS, self.testRun, accuracyPerc)

		neuronToPlot = 1
		colors = ['r']*1+['g']*1+['b']*1+['y']*1
		colors = ['blue', 'green', 'magenta', 'cyan']
		subplot(211)
		plot(M.t/ms, M.i, '.')
		legend(['A','B','C','D'], loc='upper left')			
		subplot(212)
		plot(UmM.t, UmM.v2.T/mV)	
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		if (showPlot==True):
			show()	

		return evaluateClassifier, precisionPerc

def main():
	run_gupta_paper = gupta_paper()
	evaluateClassifier, precisionPerc = run_gupta_paper.run_model()
	if evaluateClassifier: print precisionPerc
	return precisionPerc

if  __name__ =='__main__':main()