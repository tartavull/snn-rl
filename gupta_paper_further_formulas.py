from furtherFormulas.architecture_further_formulas import *
from furtherFormulas.cofactorCalculations import *
from furtherFormulas.timeAndRefracCalcs import *
from furtherFormulas.outputPrinting import *
from furtherFormulas.testingProcesses import *
from furtherFormulas.lateralInhibition import *

timeAndRefrac = timeAndRefrac	

print "Parameters that can be used when running the program are: \n\
		<trainOrTest>: if training or testing should be run\n\
		<randomization>: the amount of randomization to initialize weights for\n\
		training with.  E.x. (.5-1.0)\n\
		<posReinf>: strength of positive reinforcement.  E.x. ...\n\
		<negReinf>: strength of negative reinforcement.  E.x. ...\n\
		"

class gupta_paper:
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

	print 'initial Weights\n',W

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
			neuronIndex = self.neuronIndex
			generalClockDt = self.generalClockDt

			def __init__(self, params, *args, **kwargs):

				'''self.evaluateClassifier = True
				self.accelerateTraining = False'''

				#self.evaluateClassifier = kwargs["evaluateClassifier"]#evaluateClassifier
				#self.accelerateTraining = kwargs["accelerateTraining"]#accelerateTraining

				for key in kwargs:
					print "another keyword arg: %s: %s" % (key, kwargs[key])
				print "showing args"
				for kwarg in sys.argv:
					keyAndVal = str(kwarg).split('=')
					if (len(keyAndVal) > 1):
						currentKWArgName = str(kwarg).split('=')[0].strip()
						currentKWArgValue = str(kwarg).split('=')[1].strip()
						if currentKWArgValue == 'True':
							currentKWArgValue = True
						elif currentKWArgValue == 'False':
							currentKWArgValue = False
						
						if currentKWArgName == "evaluateClassifier":
							self.evaluateClassifier = currentKWArgValue
							# By default accelerate
							if currentKWArgValue == False: self.accelerateTraining = True
						elif currentKWArgName == "accelerateTraining":
							self.accelerateTraining = accelerateTraining

					'''b = str(kwarg).split('=')
					print len(b)
					print b[1]
					print str(kwarg).split('=')
					print str(kwarg).split('=')[0]
					#print str(kwarg).split('=')[1]'''

				NeuronGroup.__init__(self, N=dictionaryLongitude, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; v2=10*mV;UmSpikeFired=1*mV;beginRefrac=1*mV;inhibitionVoltage=prelimV',refractory=8*ms,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					neuronIndex = self.neuronIndex
					timeAndRefrac.spikeIntervalCounter = (floor(timeAndRefrac.time/timeAndRefrac.spikeIntervalUnformatted) * timeAndRefrac.spikeIntervalUnformatted)*10

					def dendCalcs(neuronIndex, ADDSObj, dendObj, spiketimes, evaluationActive):
						# Below sequentially Dirac, Tau, then Resistance are calculated every end of a spike-time interval.
						# The resulting Dend I is added to the Um calc for the ADDS soma.
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
						# dend then somaDirect calcs are done which are then used to set lat inhib.
						# Soma Um calcs are done automatically using equations entered for brian
						# once dend and somaDirect are updated
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

						'''roundedSecondsTime = math.floor(Decimal(format((ADDSObj.t), '.1f'))/Decimal(format((1.0*second), '.1f')))
						if printWeights < roundedSecondsTime:
							self.printWeights = roundedSecondsTime; printWeights(dendObj);'''
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
		ADDS = ADDSNeuronModel(self, sys.argv)			
		M = SpikeMonitor(ADDS)
		Mv = StateMonitor(ADDS, 'V', record=True)
		UmM = StateMonitor(ADDS, 'v2', record=True)
		self.M = M # for ipython compatibility
		self.UmM = UmM 
		self.weightMonitors = weightMonitors

		neuralnet.add(ADDS)
		neuralnet.add(M)
		neuralnet.add(UmM)
		self.runTime *= self.runTimeScaling # scaling factor
		neuralnet.run(self.runTime,report='text')

		OutputEvaluationResults(dend, self.testRun)

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
		show()	

	def __init__(self, *args, **kwargs):
		'''self.neuralnet = Network()
		self.dictionary = dictionary()
		self.spiketimes = dictionary.spikeTimes(self.dictionary, dictionaryLongitude, spikeInterval, testingSpikesPerChar, testingEpochs)
		self.trainingSpikeTimes = dictionary.spikeTimes(self.dictionary, dictionaryLongitude, spikeInterval, trainingSpikesPerChar, trainingEpochs)
		self.LIK = SpikeGeneratorGroup(N=15, indices=self.spiketimes[:,0], times=self.spiketimes[:,1]*ms)
		# W = W and other lines below are to avoid an odd 'reference before assignment' error
		self.W = W
		self.tauM = tauM
		self.neuronIndex = neuronIndex
		self.generalClockDt = generalClockDt
		self.runTime = runTime
		self.runTimeScaling = runTimeScaling
		self.evaluateClassifier = evaluateClassifier
		self.accelerateTraining = accelerateTraining
		self.diracScaling = diracScaling
		self.somaDirectScaling = somaDirectScaling
		self.negativeWeightReinforcement = negativeWeightReinforcement
		self.positiveWeightReinforcement = positiveWeightReinforcement
		self.timeAndRefrac = timeAndRefrac	
		self.testRun = testRun
		self.latInhibSettings = latInhibSettings'''

		'''neuralnet = Network()
		dictionary = self.dictionary()
		spiketimes = dictionary.spikeTimes(self.dictionary, dictionaryLongitude, spikeInterval, testingSpikesPerChar, testingEpochs)
		trainingSpikeTimes = dictionary.spikeTimes(self.dictionary, dictionaryLongitude, spikeInterval, trainingSpikesPerChar, trainingEpochs)
		LIK = SpikeGeneratorGroup(N=15, indices=self.spiketimes[:,0], times=self.spiketimes[:,1]*ms)
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
		latInhibSettings = latInhibSettings'''

		self.run_model()

def main(*args, **kwargs):

	if len(sys.argv) > 1:
		print "Parameters found"
		for key in kwargs:
			print "another keyword arg22: %s: %s" % (key, kwargs[key])
		print "Number of arguments111: ", len(kwargs.items())
	else:
		print "Parameters that can be used when running the program are: \n\
		<trainOrTest>: if training or testing should be run\n\
		<randomization>: the amount of randomization to initialize weights for\n\
		training with.  E.x. (.5-1.0)\n\
		<posReinf>: strength of positive reinforcement.  E.x. ...\n\
		<negReinf>: strength of negative reinforcement.  E.x. ...\n\
		"

	print "This is the name of the script: ", sys.argv[0]
	print "Number of arguments: ", len(sys.argv)
	print "The arguments are: " , str(sys.argv[1]).split('=', 1)[0]
	print "The arguments are: " , str(sys.argv[2]).split('=', 1)[0]

	run_gupta_paper = gupta_paper(sys.argv)

#if  __name__ =='__main__':main(evaluateClassifier = evaluateClassifier, negativeWeightReinforcement = negativeWeightReinforcement, positiveWeightReinforcement = positiveWeightReinforcement)
if  __name__ =='__main__':main(sys.argv)