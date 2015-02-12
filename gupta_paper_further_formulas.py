from architecture_further_formulas import *
from furtherFormulas.cofactorCalculations import *
from furtherFormulas.timeAndRefracCalcs import *
from furtherFormulas.outputPrinting import *
from furtherFormulas.testingProcesses import *
from furtherFormulas.lateralInhibition import *

timeAndRefrac = timeAndRefrac	

class gupta_paper:

	neuralnet = Network()
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	trainingSpikeTimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, 1, 300)
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	# tauD and R lines below are to avoid an odd 'reference before assignment' error
	tauD = tauD
	W = W
	R = R
	tauM = tauM
	neuronIndex = neuronIndex
	generalClockDt = generalClockDt
	runTime = runTime
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
		#dv/dt = (-v+((Rm/mV)*(SynI+DendI)))/(tauM) : volt (unless refractory)
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
			lateralInhibActive : boolean
			beginRefrac : volt
			uMScaling : volt
			preLatInhibV : volt
			inhibitionVoltage : volt
		    ''')			

		# Units are removed in dv equ below and it is converted to mV.  .001 is conversion factor
		dendriteEqs = Equations('''
			dv/dt = (((-v/mV)+((r/mV)*(w/volt)*(dirac/volt)))/(tau))*mV : volt
			V : volt
	        r : volt
	        w : volt
	        dirac : volt
	        tau : second
	        v2: volt
	        vTest : volt
	        vTest2 : volt
	        vTest3 : volt
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
			#refractoryPointCounter = timeAndRefrac.refractoryPointCounter
			calcDirac2 = None
			generalClockDt = self.generalClockDt
			def __init__(self, params):
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; v2=10*mV;UmSpikeFired=1*mV;beginRefrac=1*mV;inhibitionVoltage=prelimV',refractory=8*ms,clock=Clock(dt=self.generalClockDt))
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
						if timeAndRefrac.time == 0.000:
							dend[neuronIndex].w = W[neuronIndex]*volt
							#for inputNeuronIndex in range(numberOfPixels):
							#	dend[neuronIndex].w[inputNeuronIndex] = W[neuronIndex][inputNeuronIndex]*volt
							#print 'W::',W

						if (evaluationActive==False and timeAndRefrac.refractoryPointSwitch==True):
							# Only change weights of neuron fired
							if ADDSObj.UmSpikeFired[neuronIndex] == 1*mV:
								# Weights
								WeightChangeCalculation(neuronIndex, spiketimes, timeAndRefrac.time, self.negativeWeightReinforcement, self.positiveWeightReinforcement, M, dendObj)	
								#print 'w',' time',timeAndRefrac.time,'\t',dendObj[0].w,'\t',dendObj[1].w,'\t',dendObj[2].w,'\t',dendObj[3].w
							# Tau
							tauDCalc(neuronIndex, dendObj, self.W)
							# Resistance
							resistanceCalc(neuronIndex, dendObj, self.R, self.tauM)

						# TODO: check do I need additional loop below?
						for indexOfDend in range(dictionaryLongitude):
							ADDSObj.DendI[indexOfDend] = sum(dendObj[indexOfDend].v[:])*dendCalcScaling

					def somaDirectCalcs(neuronIndex, ADDSObj, somaDirectObj, dendObj):
						# somaDirect calcs
						# Reuse existing dendtrite values below instead of recalculating
						dotProductWandDirac =  sum(w*d for w,d in zip(dendObj[neuronIndex].w[:], dendObj[neuronIndex].dirac[:]))
						somaDirectObj.summedWandDirac[neuronIndex] = ((dotProductWandDirac*volt)/(volt**2))*self.somaDirectScaling
						#print 'test::_::somaDirectObject.summedWandDirac',somaDirectObject.summedWandDirac
						#print 'test::_::somaDirectObject.v',somaDirectObject.v
						#print 'test::_::testSomaDirect.v',testSomaDirect.v

						for neuronNumber in range(dictionaryLongitude):
							ADDSObj.SynI[neuronNumber] = somaDirectObj.v[neuronNumber]												

					def mainSimulationCalcs(ADDSObj, dendObj, somaDirectObj, spiketimes, evaluationActive):
						# dend then somaDirect calcs are done which are then used to set lat inhib.
						# Soma Um calcs are done automatically using equations entered for brian
						# once dend and somaDirect are updated
						#timeAndRefrac = self.timeAndRefrac
						preTNorm = self.timeAndRefrac.time
						tNorm = preTNorm - (floor((preTNorm/.001)*.01) * .1)
						#print 'tNorm',tNorm
						
						self.timeAndRefrac = timePeriodAndRefractoryCalcs(self.timeAndRefrac)

						#print 'time',timeAndRefrac.time,'brian time',self.t

						if (evaluationActive==True):
							intitializeTrainedModelParameters(dendObj)

						# Option to accelerate computations for training
						if self.accelerateTraining == False or (evaluationActive == False and tNorm <= .005 or tNorm >= .096):													
							if self.accelerateTraining == True and (tNorm >= .096 and tNorm < .099):
								for i in range(dictionaryLongitude):
									ADDSObj.DendI[i]=0*mV
									ADDSObj.SynI[i]=0*mV
									ADDSObj.prelimV[i]=0*mV
									ADDSObj.v[i]=0*mV#-65*mV#0*mV									
									for i2 in range(len(dend[0].v)):
										dendObj[i].v[i2] = 0*mV
									somaDirectObj.v[i] = 0*mV		

							#print 'running'
							for neuronIndex in range(dictionaryLongitude):
								dendCalcs(neuronIndex, ADDSObj, dendObj, spiketimes, evaluationActive)

								somaDirectCalcs(neuronIndex, ADDSObj, somaDirectObj, dendObj)
								# Lat inhib.
								#S.w[:]=(-1*(ADDSObj.Rm/mV*(ADDSObj.DendI[neuronIndex]+ADDSObj.SynI[neuronIndex])))*self.lateralInhibitScaling2

								#checkForResets(neuronIndex, dendObj, somaDirectObj)

							ADDSObj.v, self.latInhibSettings = lateralInhibition(ADDSObj, self.timeAndRefrac, self.latInhibSettings)

							for neuronIndex in range(dictionaryLongitude): 
								ADDSObj.v2, somaDirectObj.v2, self.timeAndRefrac = checkForResets(neuronIndex, ADDSObj, dendObj, somaDirectObj, self.timeAndRefrac)

							ADDSObj.UmSpikeFired, self.testRun = evaluateClassifierPerf(ADDSObj, self.testRun)

							'''if tNorm <= 0.003 or tNorm >= .099:
								printEvalOutputForTesting(neuronIndex, timeAndRefrac, ADDS, dend)
							else:
								sortVersionOutForTesting(neuronIndex, timeAndRefrac, ADDS, dend)'''

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
				NeuronGroup.__init__(self, N=4, model=directToSomaEqs,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					placeHolderForLaterContent = True
						#print 'DtSoma','testTime',self.t,'\tneuronIndex',neuronIndex,'\tself.v[neuronIndex]',self.v[neuronIndex],'\tsummedWandDirac',self.summedWandDirac
				self.contained_objects.append(additionToNetwork)		

		#timeAndRefrac = timeAndRefrac()

		dend = [None]*dictionaryLongitude
		testDend = [None]*dictionaryLongitude
		for firstLayerIndex in range(dictionaryLongitude):
			dend[firstLayerIndex] = DendriteNeuronModel(15)
			testDend[firstLayerIndex] = DendriteNeuronModel(15)		
		somaDirect = SomaDirectNeuronModel(4)
		testSomaDirect = SomaDirectNeuronModel(4)

		'''if self.evaluateClassifier == False:
			ADDS = ADDSNeuronModel(self)
		else:
			testADDS = ADDSNeuronModel(self)'''

		ADDS = ADDSNeuronModel(self)			

		# Synapses for lateral inhibition
		#S=Synapses(testADDS,testADDS,model='w:volt',pre='v+=w')
		'''S=Synapses(ADDS,ADDS,model='w:volt',pre='v+=w',clock=Clock(dt=self.generalClockDt))
		S.connect('i != j')
		S.w[:]=0.0*volt	'''

		M = SpikeMonitor(ADDS)
		self.M = M # for ipython compatibility
		#testM = SpikeMonitor(testADDS)
		#self.testM = testM
		Mv = StateMonitor(ADDS, 'V', record=True)
		#testMv = StateMonitor(testADDS, 'v', record=True)
		MDendI = StateMonitor(ADDS, 'DendI', record=True)
		MSynI = StateMonitor(ADDS, 'SynI', record=True)
		UmM3 = StateMonitor(ADDS, 'v2', record=True)
		self.UmM3 = UmM3 # for ipython compatibility
		#testUmM3 = StateMonitor(testADDS, 'v2', record=True)
		#self.testUmM3 = testUmM3
		somaDirectM = SpikeMonitor(somaDirect)
		somaDirectMv = StateMonitor(somaDirect, 'V', record=True)
		#synM = StateMonitor(S,'w',record=S['i!=j'])  # all synapses excluding autapses
		dendSpikeM = [None]*dictionaryLongitude
		dendVoltageM = [None]*dictionaryLongitude
		'''testDendSpikeM = [None]*dictionaryLongitude
		testDendVoltageM = [None]*dictionaryLongitude'''
		for firstLayerIndex in range(dictionaryLongitude):
			dendSpikeM[firstLayerIndex] = SpikeMonitor(dend[firstLayerIndex])
			dendVoltageM[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'v', record=True)
			'''testDendSpikeM[firstLayerIndex] = SpikeMonitor(dend[firstLayerIndex])
			testDendVoltageM[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'v', record=True)'''

			neuralnet.add(dend[firstLayerIndex])
			neuralnet.add(dendSpikeM[firstLayerIndex])
			neuralnet.add(dendVoltageM[firstLayerIndex])
			'''neuralnet.add(testDend[firstLayerIndex])
			neuralnet.add(testDendSpikeM[firstLayerIndex])
			neuralnet.add(testDendVoltageM[firstLayerIndex])	'''		

		neuralnet.add(ADDS)
		#neuralnet.add(testADDS)
		neuralnet.add(M)
		#neuralnet.add(testM)
		#neuralnet.add(testMv)
		neuralnet.add(UmM3)
		#neuralnet.add(testUmM3)
		neuralnet.add(somaDirect)
		#neuralnet.add(testSomaDirect)
		neuralnet.add(somaDirectM)
		neuralnet.add(somaDirectMv)
		#neuralnet.add(S)

		#neuralnet.add(synM)
		#neuralnet.run(70*ms,report='text')
		self.runTime *= 10 # scaling factor
		neuralnet.run(self.runTime,report='text')

		OutputEvaluationResults(dend, self.testRun)

		neuronToPlot = 1
		colors = ['r']*1+['g']*1+['b']*1+['y']*1
		colors = ['blue', 'green', 'magenta', 'cyan']
		subplot(221)
		plot(M.t/ms, M.i, '.')
		#plot(M.t/ms, M.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(222)
		#plot(testM.t/ms, testM.i, '.')
		#plot(testM.t/ms, testM.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(223)
		#plot(UmM3.t, UmM3.v2.T/mV)
		#plot(UmM3.t, UmM3.v2.T/mV)
		#plot(testMv.t, testMv.v.T/mV)	
		plot(Mv.t, Mv.V.T/mV)	
		legend(['A','B','C','D'], loc='upper left')			
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		subplot(224)
		#plot(testUmM3.t, testUmM3.v2.T/mV)	
		plot(UmM3.t, UmM3.v2.T/mV)	
		#plot(testUmM3.t, testUmM3.v2Test.T/mV)
		#plot(testMv.t, testMv.vTest.T/mV)
		#plot(dendVoltageM[0].t, dendVoltageM[0].v.T/mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		'''subplot(225)
		plot(somaDirectM.t/ms, somaDirectM.i, '.')
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')'''
		'''subplot(226)
		plot(somaDirectMv.t, somaDirectMv.v2.T/mV) 	
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')'''
		'''subplot(225)
		plot(synM.t / ms, synM[0, :].w / nS) 	
		#plot(testUmM3.t, testUmM3.v2Test.T/mV)
		#plot(testMv.t, testMv.vTest.T/mV)
		#plot(dendVoltageM[0].t, dendVoltageM[0].v.T/mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')'''
		show()	

	def __init__(self):
		self.run_model()

def main():
	run_gupta_paper = gupta_paper()

if  __name__ =='__main__':main()
