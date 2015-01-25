from architecture_further_formulas import *

class gupta_paper:
	neuralnet = Network()
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	testSpiketimes = spiketimes
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	neuronIndex = 0
	time = 0.0000
	testTime = 0.0000
	SummedDendriteGroup = 0;
	SynapseToSoma = 0;
	SigmaWDSyn = .0375 * mV
	SigmaWDSynMod = .5
	RmEq = Rm #Rm * mV
	DendR = .0375 * mV
	DendW = 1
	DendDirac = 1	
	# tauD and R lines below are to avoid an odd 'reference before assignment' error
	tauD = tauD
	testTauD = testTauD
	R = R
	testR = testR
	tauM = tauM
	lastSpikeInterval = -0.001#0.0
	testLastSpikeInterval = -0.001#0.0
	currentEpoch = 0
	timeStepInterval = 0.001
	refractoryPointCounter = 0.0
	testRefractoryPointCounter = 0.0
	refractoryPointSwitch = False
	testRefractoryPointSwitch = False
	neuronIndexCounter = 0.0
	neuronIndexSwitch = False
	#timing = []
	truePositiveSpikeResults = 0
	falsePositiveSpikeResults = 0
	trueNegativeSpikeResults = 0
	falseNegativeSpikeResults = 0
	totalSpikeIntervals = 12	
	testSpikesFiredInInterval = np.array([[False]*(totalSpikeIntervals+1)]*dictionaryLongitude)
	IdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	testIdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	IsSpikeFired = np.array([False]*dictionaryLongitude); 
	testIsSpikeFired = np.array([False]*dictionaryLongitude);	
	UmSpikeFired = np.array([False]*dictionaryLongitude); 
	#testUmSpikeFired = np.array([False]*dictionaryLongitude); 
	IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	testIdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	IsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	testIsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	UmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	testUmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	refractoryGraphingSwitch = False

	correctSpikes = np.array([[1]*totalSpikeIntervals]*dictionaryLongitude);#[[1] * totalSpikeIntervals]*dictionaryLongitude
	correctSpikes[0][3:12] = 0
	correctSpikes[1][0:3] = 0
	correctSpikes[1][6:12] = 0
	correctSpikes[2][0:6] = 0
	correctSpikes[2][9:12] = 0
	correctSpikes[3][0:9] = 0
    
	M = None
	testM = None
	UmM3 = None
	testUmM3 = None
	testSpikeIntervalCounter = 0.0

	print 'initial Weights\n',W

	def run_model(self):
		neuralnet = self.neuralnet
		dictionary = self.dictionary
		spiketimes = self.spiketimes
		#self.ADDSNeuronModel2 = self.ADDSNeuronModel
		#dv/dt = (-v+((RmZ/volt)*(IdVZ)))/(tauMZ) : volt (unless refractory)
		#dv/dt = (-v+((RmZ/volt)*(IdVZ)))/(tauMZ) : volt (unless refractory)
		eqs = Equations('''
			dv/dt = (-v+((RmZ/mV)*(IdVZ)))/(20*ms) : volt (unless refractory)
			RmZ = 80*mV : volt
			tauMZ = 30*ms : second
			IdVZ : volt
	        V : volt
	        DendI : volt
	        SynI : volt
	        v2 : volt	
			testUmSpikeFired : volt	
			lateralInhibActive : boolean
			refractoryPointSwitch : boolean	
		    ''')			

		#dv/dt = (-v+((r/volt)*dirac))/(tau) : volt (unless refractory)
		#dv/dt = (-v+((r/mV)*(w/mV)*dirac))/(tau) : volt (unless refractory)
		#dv/dt = (-v+(3*mV))/(20*ms) : volt (unless refractory)
		#dv/dt = (-v+((r/volt)*(w/volt)*dirac*.001))/(tau) : volt (unless refractory)
		dendriteEqs = Equations('''
			dv/dt = (-v+((r/volt)*(w/volt)*dirac))/(tau) : volt (unless refractory)
			V : volt
	        r : volt
	        w : volt
	        dirac : volt
	        tau : second
	        v2: volt
			''')

		directToSomaEqs = Equations('''
			dv/dt = (-v+summedWandIDirac)/(tauMZ) : volt (unless refractory)
			tauMZ = 30*ms : second
			IdVZ : volt
	        V : volt
	        DendI : volt
	        SynI : volt
	        v2 : volt
			''')		

		'''myclock=Clock(dt=1*ms)
		testMyclock=Clock(dt=1*ms)		

		ADDS = NeuronGroup(N=4, model=eqs,threshold='v>.01*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV',refractory=10*ms,clock=myclock)
		testADDS = NeuronGroup(N=4, model=testEqs,threshold='vTest>.01*mV', reset='vTest=-0.002 * mV; dvTest=0; IdVZTest = 0*mV;v2Test=10*mV;lateralInhibActive=True',refractory=10*ms,clock=testMyclock)'''

		class ADDSNeuronModel(NeuronGroup, gupta_paper): 
			neuronIndex = self.neuronIndex
			refractoryPointCounter = self.refractoryPointCounter
			calcDirac2 = None
			def __init__(self, params): 
				myclock=Clock(dt=1*ms)
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV',refractory=10*ms,clock=myclock)
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.002*mV', reset='v=-0.002 * mV',refractory=0.3*ms) 
				@network_operation 
				def additionToNetwork(): 
					neuronIndex = self.neuronIndex
					spikeIntervalCounter = (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10
					#print 'neuronIndex', self.neuronIndex, 'self.time', self.time,'ADDS.t', ADDS.t, 'ADDS.v', ADDS.v[neuronIndex]
					#print 'neuronIndex', neuronIndex, 'self.time', self.time, 'refractoryPointCounter',self.refractoryPointCounter,'ADDS.t', ADDS.t, 'testADDS.vTest', testADDS.vTest
					#print 'neuronIndex', self.neuronIndex, 'self.time', self.time, 'ADDS.t', ADDS.t, 'ADDS.IdVZ', ADDS.IdVZ[neuronIndex]
					#print 'neuronIndex', neuronIndex, 'self.time', self.time, 'ADDS.t', ADDS.t, 'testADDS.IdVTest', testADDS.IdVTest

					#returnUm(self)

					def timePeriodAndRefractoryCalcs():
						neuronIndex = self.neuronIndex
						refractoryPointCounter = self.refractoryPointCounter
						#refractoryPointSwitch = self.refractoryPointSwitch
						neuronIndexCounter = self.neuronIndexCounter
						neuronIndexSwitch = self.neuronIndexSwitch
						time = self.time

						time = time + self.timeStepInterval
						#self.timing.append(time)
						refractoryPointCounter = refractoryPointCounter + self.timeStepInterval
						# if statement below corrects for offset of spiketimes starting at .1 sec.
						if time >= .1:
							neuronIndexCounter = neuronIndexCounter + self.timeStepInterval

						#print '===+++__','time',time,'refractoryPointCounter',refractoryPointCounter

						# At the end of each spike time interval refractory period is turned off and weight changes
						# are calculated.  Refractory turning off here * is not correct * because it can occur for
						# less than the set refractory period.  I am just using it for the time being for testing.
						self.refractoryPointSwitch = False
						if refractoryPointCounter >= spikeIntervalUnformatted:
							#refractoryPointCounter = 0.009
							refractoryPointCounter = 0.000
							self.refractoryPointSwitch = True
							self.lastSpikeInterval = time-self.timeStepInterval
							#self.lastSpikeInterval = self.time
							WeightChangeCalculation()
							# Report weights at different times
							#if self.time == 1.081 or self.time == 2.081 or self.time == 4.081 or self.time == 8.081 or self.time == 12.081:
							'''if self.currentEpoch == 2 or self.currentEpoch == 4 or self.currentEpoch == 8 or self.currentEpoch == 10 or self.currentEpoch == 12:
								print 'self.time:\t',self.time,'\tself.currentEpoch\t',self.currentEpoch,'\tcurrent neuront\t',self.neuronIndex
								print '0',W[0][:]
								print '1',W[1][:]
								print '2',W[2][:]
								print '3',W[3][:]'''

						neuronIndexSwitch = False
						if neuronIndexCounter >= (spikeIntervalUnformatted*3):
							#neuronIndexCounter = 0.009
							neuronIndexCounter = 0.0000
							neuronIndexSwitch = True
							# add epoch every epoch ms duration.  Occurs when a full series of input charactors has
							# occured and a new one starts.  self.neuronIndex == 3 really means it will roll back to
							# self.neuronIndex == 3 after neuronIndex is added
							if self.neuronIndex == 3:
								self.currentEpoch = self.currentEpoch + 1		
							# Clear prior settings.  TODO: figure out if this would be updated anyway if dirac is triggered by spiketimes
							# if the spike input would have really activated it.
							dend[neuronIndex].dirac = [0*volt]*len(dend[neuronIndex][:])

						# revision of when new neuron index is used.  new neuron every 300ms
						if neuronIndexSwitch == True:
							neuronIndex = neuronIndex + 1

						if neuronIndex == 4: 
							neuronIndex = 0

						for vCheck in range(len(ADDS.v2)):
							if ADDS.v2[vCheck] == 10*mV:
								ADDS.v2[vCheck] = 0*mV
								# Clear prior settings.  These are workarounds to avoid an unwanted lingering voltage in
								# the dend input after the training section for its neuron index is done.  TODO: see if there
								# is a better biologically accurate way for this if needed
								dend[neuronIndex].dirac = [0*volt]*len(dend[neuronIndex][:])
								dend[neuronIndex].v = [0*volt]*len(dend[neuronIndex][:])

						self.neuronIndex = neuronIndex
						self.time = time
						self.refractoryPointCounter = refractoryPointCounter
						#self.refractoryPointSwitch = refractoryPointSwitch
						self.neuronIndexCounter = neuronIndexCounter			
						self.neuronIndexSwitch = neuronIndexSwitch

					def WeightChangeCalculation():
						for inputNeuronIndex in range(numberOfPixels):
							## General STDP learning rule implementation ##
							# Find SpikePreSyn if it exists
							# Note: could replace this with use of dictionary data structure for lookups if convenient
							# for processing time later

							# SpikePreSyn not found than make it max distance from SpikePostSyn
							# TODO SpikePreSyn = 0.1 # (Max spike time interval distance)
							# round up to nearest 100ms interval
							SpikePreSyn = (math.ceil(self.time*10)*.1)*second
							# SpikePostSyn not found than make it max distance from SpikePreSyn
							# TODO SpikePostSyn = 0.0
							#SpikePostSyn = 0*ms
							SpikePostSyn = (math.floor(self.time*10)*.1)*second

							preSynSpikeFound = False
							postSynSpikeFound = False

							spikeCollection = spiketimes
							NumberOfSpikes = shape(spikeCollection)[0]
							for i in range(NumberOfSpikes):
								CurrentSpikeNueron = spikeCollection[i][0]
								CurrentSpikeTime = spikeCollection[i][1]*second

								# exit loop once values below current time elapsed have all been checked
								if CurrentSpikeTime > (self.time*second):
									#print 'eeww1:\tspiketimes\t',spiketimes,'\tqq\t',CurrentSpikeTime
									break

								# (self.time*second-.1*second) is to check if in relevant time window below.
								# Note: that may not be a good cut off and I should check it
								# added < self.time just for testing
								if CurrentSpikeNueron == inputNeuronIndex and CurrentSpikeTime > (self.time*second-.1*second):# and CurrentSpikeTime < (self.time*second):
									SpikePreSyn = CurrentSpikeTime
									preSynSpikeFound = True

							# Find SpikePostSyn if exists
							#spikeCollection = M.spikes
							spikeCollection = M.it
							#print 'spikeCollection',spikeCollection
							#print 'spikeCollection', spikeCollection
							NumberOfSpikes = len(spikeCollection[0])
							for i in range(NumberOfSpikes):
								#CurrentSpikeNueron = spikeCollection[i][0]
								#CurrentSpikeTime = spikeCollection[i][1]*1000
								CurrentSpikeNueron = spikeCollection[0][i]
								#CurrentSpikeTime = spikeCollection[1][i]*1000
								CurrentSpikeTime = spikeCollection[1][i]*10
								#print 'CurrentSpikeNueron',CurrentSpikeNueron,'CurrentSpikeTime',CurrentSpikeTime,'self.time*second',self.time*second

								# exit loop once values below current time elapsed have all been checked
								# Disabled due to spikeCollection not being sorted and causing break too early
								#if CurrentSpikeTime > (self.time*second):
								#	break

								# (self.time*second-.1*second) is to check if in relevant time window below.
								# Note: that may not be a good cut off and I should check it
								# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
								# added (.1+.008)
								if CurrentSpikeNueron == self.neuronIndex and CurrentSpikeTime >= (self.time*second-(.1+.008)*second) and CurrentSpikeTime <= (self.time*second):
									SpikePostSyn = CurrentSpikeTime	
									postSynSpikeFound = True	

							if preSynSpikeFound == False or postSynSpikeFound == False:
								# Todo: watch more for spikes in a relevant time frame but are actually from different synapses and would
								# give different results?  could cause unwanted results? brain would know synapse difference?

								## todo: currently included is logic for processing if SpikePostSyn or SpikePreSyn were not found?

								if preSynSpikeFound == False and postSynSpikeFound == True:
									SpikePreSyn = SpikePostSyn + .1*second
									#print 'F T'
								elif preSynSpikeFound == True and postSynSpikeFound == False:
									# Note: does postSynSpikeFound need to be true for the weight to be allowed to increase?
									SpikePostSyn = SpikePreSyn - .1*second
									#print 'T F'
								#else:
									#print 'F F'

							#SpikePreSyn = 500*ms
							#SpikePostSyn = 700*ms
							#print 'neuronIndex', self.neuronIndex, 'inputNeuronIndex',inputNeuronIndex,'self.time', self.time, 'ADDS.v2',ADDS.v2, 'preSynSpikeFound',preSynSpikeFound,'postSynSpikeFound',postSynSpikeFound,'SpikePostSyn', SpikePostSyn, 'SpikePreSyn', SpikePreSyn

							# Find DeltaW
							DeltaW = returnDeltaW(SpikePreSyn, SpikePostSyn)  

							# Find new weight
							WOld = W[self.neuronIndex][inputNeuronIndex];

							W[self.neuronIndex][inputNeuronIndex] = returnNewW(WOld, DeltaW);

							dend[self.neuronIndex].w[inputNeuronIndex] = W[self.neuronIndex][inputNeuronIndex]*volt

					def returnDeltaW(SpikePreSyn, SpikePostSyn):
						DeltaSpikeTime = SpikePreSyn - SpikePostSyn
						# TODO: figure out if DeltaW = 0 is really fine for init value
						DeltaW = 0
						# changed below line from that content in the article for a workaround/temp fix
						#if DeltaSpikeTime < 0:
						if DeltaSpikeTime <= 0:
							DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus*ms)))
							# for testing
							DeltaW = 1
						elif DeltaSpikeTime > 0:
							DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus*ms)))	
							# for testing
							DeltaW = -1

						# testing
						#print 'self.time',self.time,'SpikePreSyn',SpikePreSyn,'SpikePostSyn',SpikePostSyn,'DeltaSpikeTime',DeltaSpikeTime,'DeltaW',DeltaW
						#if 	(t*ms) - SpikePreSyn > (1*ms):
						#	DeltaW = .5
						#else:
						#	DeltaW = -.5				

						return DeltaW

					def returnNewW(WOld, DeltaW):
						# W is weight
						synapseType = 'excitatory'
						# Default 'excitatory' values
						WMin = 0;
						WMax = 1;
						if synapseType == 'inhibitory':
							WMin = -1;
							WMax = 0;		
							
						#print 'DeltaW', DeltaW
						## Relaxation rule implementation ##
						if DeltaW >= 0:				
							WNew = WOld + (LearningRate * (DeltaW * (WMax - WOld)))
							
						elif DeltaW < 0:	
							WNew = WOld + (LearningRate * (DeltaW * (WOld - WMin)))
						#print 'WOld', WOld
						#print 'WNew', WNew

						return WNew;
						##	

					def tauDCalc(neuronIndex, dendObject, W):
						# Weight loop
						for WIndex in range(len(W[0][:])):
							if abs(W[neuronIndex][WIndex]) <= 1:
								dendObject[neuronIndex].tau[WIndex] = (tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)) * .001 * second# .001 is scaling factor found in dirac method

						return dendObject[neuronIndex].tau

					def resistanceCalc(neuronIndex, dendObject, R):
						tauM = self.tauM
						#Resistance loop
						for RIndex in range(len(R[0][:])):				
							if (dendObject[neuronIndex].tau[RIndex]*.001 == tauM*second):
								# avoid a division by 0 issue
								tauM = tauM - .000001

							dendObject[neuronIndex].r[RIndex] = (((((dendObject[neuronIndex].tau[RIndex]*.001)/ms)*neuronFiringThreshold) / Rm) * ((tauM / ((dendObject[neuronIndex].tau[RIndex]/ms)*.001) ) ** (tauM / (tauM - ((dendObject[neuronIndex].tau[RIndex]/ms)*.001) ))))*volt

						return dendObject[neuronIndex].r				

					def mainSimulationCalcs():
						# dirac
						#dend[neuronIndex].dirac = self.diracCalc(dend, False)#super(gupta_paper, run_model).diracCalc(dend, False)
						dend[neuronIndex].dirac = self.diracCalc(dend, False, self.time, math.e, self.neuronIndex, spiketimes,self.lastSpikeInterval)
						'''print 'neuronIndex',self.neuronIndex,'time',self.time
						print 'dend.dirac',dend[0].dirac,dend[1].dirac,dend[2].dirac,dend[3].dirac
						print 'dend.v',dend[0].v,dend[1].v,dend[2].v,dend[3].v
						print 'sum(dend[neuronIndex].v[:])',sum(dend[neuronIndex].v[:])
						print 'Prelim Weights\n0',dend[0].w[:]
						print '1',dend[1].w[:]
						print '2',dend[2].w[:]
						print '3',dend[3].w[:]
						print 'Prelim Tau\n0',dend[0].tau[:]
						print '1',dend[1].tau[:]
						print '2',dend[2].tau[:]
						print '3',dend[3].tau[:]					
						print 'Prelim Res\n0',dend[0].r[:]
						print '1',dend[1].r[:]
						print '2',dend[2].r[:]
						print '3',dend[3].r[:]'''
						'''for i in range(4):
							if i != neuronIndex:
								dend[i].dirac = [0*volt]*15'''

						# Calculate tauD
						tauD = self.tauD
						tauD = tauDCalc(neuronIndex, dend, W)
						#self.tauD = tauD # this line is just for resutls reporting and can be removed later
						#print 'tauD',tauD

						# Calculate resistance
						R = self.R
						R = resistanceCalc(neuronIndex, dend, R)
						#self.R = R # this line is just for resutls reporting and can be removed later
						#print 'R',R

						# Record sum of dend voltage
						#ADDS.IdVZ[neuronIndex] = sum(dend[neuronIndex].v[:])#0.0375*mV 
						for i in range(4):
							ADDS.IdVZ[i] = sum(dend[i].v[:])
						'''for i in range(4):
							if i != neuronIndex:
								ADDS.IdVZ[i] = 0*volt'''
						#print '::ADDS.IdVZ::\t', ADDS.IdVZ
						#print 'ADDS.v2',ADDS.v2

						timePeriodAndRefractoryCalcs()

					mainSimulationCalcs()

				self.contained_objects.append(additionToNetwork)	

			def diracCalc(self, IDend, evaluationActive, time, e, neuronIndex, spiketimes,lastSpikeInterval):
				# This method only calculates dirac
				'''time = self.time
				e = math.e
				neuronIndex = self.neuronIndex'''


				#print ':::range(len(IDend[neuronIndex][:])\t','self.time\t',self.time,'\tneuronIndex\t',neuronIndex,'\t',IDend,'\t',IDend[neuronIndex]

				# normalize t to count just time within spike interval.  Formula from article just uses t in a way
				# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
				#tNorm = time - (floor(time/spikeIntervalUnformatted) * spikeIntervalUnformatted) 

				dendGroup = [None]*len(dend[neuronIndex][:])
				for IdIndex in range(len(dend[neuronIndex][:])):
					#tauDen = tauD[neuronIndex][IdIndex]
					#r = R[neuronIndex][IdIndex]
					#w = W[neuronIndex][IdIndex]

					# set tpresyn initially too far to trigger dirac
					#tPreSyn = -t - 1
					tPreSyn = time + 1
					for presynInd in range(shape(spiketimes)[0]):
						comparedSpikeTime = spiketimes[presynInd][1]

						#print '+++===self.time',self.time,'neuronIndex',neuronIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0],'comparedSpikeTime',comparedSpikeTime,'IdIndex',IdIndex,'lastSpikeInterval',lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted#,'self.lastSpikeInterval',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'IdIndex',IdIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0]
						# this is true? -> commented out logic below due to spiketimes perhaps not being in sorted order, a relevant entry could occur after an irrelevant one
						#if (comparedSpikeTime-spikeIntervalUnformatted) > t:
						# spike times are in order of earliest to latest times grouped with input pixel number.  therefore the below cutoff in the loop is fine.
						#if comparedSpikeTime > self.lastSpikeInterval:
						if comparedSpikeTime > (lastSpikeInterval + spikeIntervalUnformatted):
							break

						#print ':::',comparedSpikeTime,self.lastSpikeInterval,spikeIntervalUnformatted
						##print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'spiketimes[presynInd][0] == IdIndex',(spiketimes[presynInd][0] == IdIndex),'IdIndex',IdIndex,'comparedSpikeTime > self.lastSpikeInterval',(comparedSpikeTime > self.lastSpikeInterval),'comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted)',(comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted))
						'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
						# 0.0 to 300 ms, as the input is intended
						# changed <= to < below to see if that works better
						if spiketimes[presynInd][0] == IdIndex and (comparedSpikeTime-spikeIntervalUnformatted) < t:'''
						# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
						#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted):
						#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - (spikeIntervalUnformatted*2)) and comparedSpikeTime <= (self.lastSpikeInterval - spikeIntervalUnformatted):
						#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted) and comparedSpikeTime <= self.lastSpikeInterval:
						if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > lastSpikeInterval and comparedSpikeTime <= (lastSpikeInterval + spikeIntervalUnformatted):
							## adding below if statement to avoid constant input stimulus in evaluation, that seem more accurate to how the article did it's evaulation.
							#print 'tNorm',tNorm
							tNorm2 = time - (floor((time/.001)*.01) * .1)
							#print 'tNorm2',tNorm2
							#if tNorm2 < .02 or evaluationActive == False:
							if True:
								tPreSyn = comparedSpikeTime
								#print '++++++++++'

					#Id2 = IDend[neuronIndex].v[IdIndex]
					#Dt = time - tPreSyn
					Dt = Decimal(format(time, '.8f')) - Decimal(format(tPreSyn, '.8f'))
					#print '++++++++Dt++++++++',Dt,'Decimal(format(time, .8f))',Decimal(format(time, '.8f')),'Decimal(format(tPreSyn, .8f))',Decimal(format(tPreSyn, '.8f'))
					# self.time == lines below are a workaround for initialization values
					if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
						Dt = 1.0
					##print 'neuronIndex',neuronIndex,'tPreSyn000:',tPreSyn,'Dt',Dt
					# dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
					# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
					# though
					# simplify dirac for testing
					#DiracFun = 0
					#DiracFun = 1/Dt

					#SpikeModCoeff = (r*w*DiracFun)

					# dirac test
					# t in dirac forumula means curent time or last spike time?		
					#if (t > -(Dt/2) and t < (Dt/2)):
					#if Dt <= spikeIntervalUnformatted:
					#SpikeModCoeff = (r*w)
					SpikeModCoeff = 0 # leaving this for reference
					if evaluationActive == True: SpikeModCoeff *= 0.018#0.018#0.015#0.003#0.01#0.0206325
					
					#print '***ni:\t',neuronIndex,'\tIdIndex\t',IdIndex,'\t***w:\t',w

					# correct for scaling
					#tauDen = tauDen * .001
					
					#if -Dt/2<(tNorm or t?)<Dt/2:
					#if Dt <= 0.0:
					if Dt >= 0.0:
						if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
							DiracFun = 1.00
						#SpikeModCoeff = .011
						#DiracFun = 0.66
						#SpikeModCoeff = (SpikeModCoeff*DiracFun)
						#tauDen = .03	
						#IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
						#if evaluationActive == True: print 'test IDend2:\tneuronIndex',neuronIndex,'IdIndex',IdIndex,'self.time',self.time,'IDend[neuronIndex]',IDend[neuronIndex],'SpikeModCoeff',SpikeModCoeff
						if Dt != 0: 
							DiracFun = 1/Dt
						else:
							DiracFun = 1000
						dendGroup[IdIndex] = float(DiracFun)*volt*.001 #add .001 as scaling factor
					else:
						DiracFun = 0
						#SpikeModCoeff = (SpikeModCoeff*DiracFun)															
						#IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
						dendGroup[IdIndex] = float(DiracFun)*volt*.001						

					#print 'neuronIndex',neuronIndex,'IdIndex',IdIndex,'dirac\t',DiracFun
				return dendGroup				

		class DendriteNeuronModel(NeuronGroup): 
			def __init__(self, params): 
				dendClock=Clock(dt=1*ms)
				NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>10*mV', reset='v=-0.002*mV; dv=0; v2=10*mV',refractory=10*ms,clock=dendClock) 
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.01*mV', reset='v=-0.002*mV; dv=0; v2=10*mV',refractory=10*ms,clock=dendClock) 
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.002*mV', reset='v=-0.002 * mV',refractory=0.3*ms) 
				@network_operation 
				def additionToNetwork(): 
					placeHolderForLaterContent = True
				self.contained_objects.append(additionToNetwork)

		dend = [None]*dictionaryLongitude
		testDend = [None]*dictionaryLongitude
		for firstLayerIndex in range(dictionaryLongitude):
			dend[firstLayerIndex] = DendriteNeuronModel(15)
			testDend[firstLayerIndex] = DendriteNeuronModel(15)					

		class testADDSNeuronModel(NeuronGroup, gupta_paper): 
			testW = None
			testR = None
			#time = self.time
			#intitializeParameters()

			def __init__(self, params): 
				#self.time = time
				testMyclock=Clock(dt=1*ms)
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV;testUmSpikeFired=1*mV;lateralInhibActive=True',refractory=10*ms,clock=testMyclock)
				#NeuronGroup.__init__(self, N=4, model=testEqs,threshold='vTest>.002*mV', reset='vTest=-0.002 * mV',refractory=0.3*ms) 
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.002*mV', reset='v=-0.002 * mV',refractory=0.3*ms) 
				@network_operation 
				def additionToNetwork(): 
					timeAndRefractoryCalcs()

					mainTestSimulationCalcs()			
					#print '*******time********',ADDS.time					
					# classifier performance test
					evaluateClassifier()				

				def mainTestSimulationCalcs():
					# for each charactor input below test the classfication results
					test1 = [None]*4
					for neuronIndex in range(dictionaryLongitude):
						#testDend[neuronIndex].dirac = ADDS.diracCalc(testDend, False)
						testDend[neuronIndex].dirac = ADDS.diracCalc(testDend, True, self.testTime, math.e, neuronIndex, spiketimes,self.testLastSpikeInterval)
						#print '^^^diracCalc^^^',False, self.testTime, math.e, neuronIndex, self.testLastSpikeInterval

						for indexOfDend in range(dictionaryLongitude):
							testADDS.IdVZ[indexOfDend] = sum(testDend[indexOfDend].v[:])

					for vCheck in range(len(testADDS.v2)):
						if self.v2[vCheck] == 10*mV:
							self.v2[vCheck] = 0*mV	

					'''print 'testDend.dirac',testDend[0].dirac,testDend[1].dirac,testDend[2].dirac,testDend[3].dirac
					print 'testDend.v',testDend[0].v,testDend[1].v,testDend[2].v,testDend[3].v
					print 'sum(testDend[neuronIndex].v[:])',sum(testDend[neuronIndex].v[:])							
					print 'Prelim Weights\n0',testDend[0].w[:]
					print '1',testDend[1].w[:]
					print '2',testDend[2].w[:]
					print '3',testDend[3].w[:]
					print 'Prelim Tau\n0',testDend[0].tau[:]
					print '1',testDend[1].tau[:]
					print '2',testDend[2].tau[:]
					print '3',testDend[3].tau[:]					
					print 'Prelim Res\n0',testDend[0].r[:]
					print '1',testDend[1].r[:]
					print '2',testDend[2].r[:]
					print '3',testDend[3].r[:]	'''				

				def timeAndRefractoryCalcs():
					self.testTime = self.testTime + self.timeStepInterval
					self.testRefractoryPointCounter = self.testRefractoryPointCounter + self.timeStepInterval	
					self.testSpikeIntervalCounter = (floor(self.testTime/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10

					#print 'self.testSpikeIntervalCounterself.testSpikeIntervalCounter=',self.testSpikeIntervalCounter,'self.testTime',self.testTime,'spikeIntervalUnformatted',spikeIntervalUnformatted

					self.testRefractoryPointSwitch = False
					#print 'testRefractoryPointCountertestRefractoryPointCounter=',self.testRefractoryPointCounter
					if self.testRefractoryPointCounter >= spikeIntervalUnformatted:
						self.testRefractoryPointCounter = 0.000
						self.testRefractoryPointSwitch = True
						self.testLastSpikeInterval = self.testTime-self.timeStepInterval	

					#print '__+++===','time',self.testTime,'refractoryPointCounter',self.testRefractoryPointCounter

				def evaluateClassifier():
					# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
					# the +1 in (self.timeotalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
					if self.testSpikeIntervalCounter < (self.totalSpikeIntervals+1):
						for neuronIndex in range(dictionaryLongitude):				
							# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
							# (self.testSpikeIntervalCounter-1) is to correct for refractoryPointSwitch occuring after spikeInterval it addresses.
							#print 'refractoryPointSwitch = ', self.refractoryPointSwitch
							#print 'refractoryPointSwitch =', self.testRefractoryPointSwitch
							#if self.refractoryPointSwitch == true and (self.testSpikeIntervalCounter > 0):
							if self.testRefractoryPointSwitch == True and (self.testSpikeIntervalCounter > 0):
								#if self.testSpikesFiredInInterval[neuronIndex][self.testSpikeIntervalCounter-1] == False:


								if (testADDS.testUmSpikeFired[neuronIndex] == 1*mV) and (self.testSpikeIntervalCounter < self.totalSpikeIntervals):
									if (self.correctSpikes[neuronIndex][self.testSpikeIntervalCounter] == 1):
										self.truePositiveSpikeResults = self.truePositiveSpikeResults + 1	
									else:
										self.falsePositiveSpikeResults = self.falsePositiveSpikeResults + 1	
									#self.testSpikesFiredInInterval[neuronIndex][self.testSpikeIntervalCounter] = True	
								elif (testADDS.testUmSpikeFired[neuronIndex] == 0*mV) and (self.testSpikeIntervalCounter < self.totalSpikeIntervals):
									if (self.correctSpikes[neuronIndex][(self.testSpikeIntervalCounter-1)] == 1):
										self.falseNegativeSpikeResults = self.falseNegativeSpikeResults + 1		
									else:
										self.trueNegativeSpikeResults = self.trueNegativeSpikeResults + 1	
								testADDS.testUmSpikeFired[neuronIndex] = 0*mV	
						#print 'results',self.truePositiveSpikeResults,self.falsePositiveSpikeResults,self.trueNegativeSpikeResults,self.falseNegativeSpikeResults
				self.contained_objects.append(additionToNetwork)

			def intitializeParameters():
				#tNorm = self.time - (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
				tNorm = self.testTime - (floor((self.testTime/.001)*.01) * .1)

				'''testW = [[1, 0, 1,0, 1, 0,0, 0, 0,0, 1, 0,0, 1, 0],
				[0, 0, 1,0, 1, 0,0, 0, 1,0, 1, 0,0, 0, 1],
				[0, 0, 0,0, 1, 1,0, 1, 1,0, 1, 1,0, 0, 0],
				[0, 0, 1,0, 1, 0,0, 1, 0,0, 1, 0,0, 0, 1]]'''
				testW = [[0.7004951, 0.00299856, 0.70085791, 0.00369599, 0.70079641, 0.00245995, 0.00433145, 0.00302166, 0.00261059, 0.00338837, 0.70107221, 0.00382794, 0.00391812, 0.7002158, 0.00261986],
				[0.00355271, 0.00352117, 0.70031808, 0.00371104, 0.70217063, 0.00440704, 0.00331328, 0.00297992, 0.70041271, 0.00235294, 0.70154393, 0.00271551, 0.00442378, 0.00403136, 0.70143686],
				[0.00401543, 0.00265892, 0.00361079, 0.00434494, 0.66817783, 0.66848986, 0.00432942, 0.66698178, 0.66919802, 0.00411504, 0.66889164, 0.66714679, 0.00295245, 0.00363973, 0.00374949],
				[0.00527166, 0.00558407, 0.70225087, 0.00466969, 0.70021103, 0.00467829, 0.00430914, 0.69986055, 0.00358859, 0.00480534, 0.70199339, 0.00578071, 0.00583865, 0.00618692, 0.70014918]]

				self.testW = testW

				# Calculate tauD
				#testTauD = self.testTauD
				#testTauD = tauDCalc(neuronIndex, testTauD, testW)
				#print 't:',self.time,'neuronIndex',neuronIndex,'main calc tauDCalc()', testTauD

				testTauD = [[10.37748461, 29.90900423, 10.33746184, 29.88206559, 10.37939543, 29.928977, 29.932435, 29.89185801, 29.87504059, 29.87123073, 10.39396579, 29.90699629, 29.87304979, 10.39251402, 29.89803348],
				[29.91250761, 29.91511768, 10.36228604, 29.92972817, 10.37676875, 29.92606665, 29.8806098, 29.93119539, 10.37107451, 29.92359325, 10.33742747, 29.91957909, 29.89789419, 29.93362677, 10.38618516],
				[29.89507952, 29.84756941, 29.89640735, 29.85177368, 12.28555811, 12.29181016, 29.88055016, 12.3308681, 12.3610359, 29.90614938, 12.35836972, 12.31503375, 29.87499997, 29.911543, 29.88966829],
				[29.89320623, 29.86032567, 10.35528394, 29.85394614, 10.32808049, 29.85050638, 29.88449768, 10.33466202, 29.82325308, 29.87764556, 10.38815977, 29.83505713, 29.90864456, 29.85475254, 10.32299042]]

				# Calculate resistance
				#testR = self.testR
				#testR = resistanceCalc(neuronIndex, testTauD, testR)

				#Sample output for testing
				'''testR = [[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463]]'''
				testR = [[0.0065761, 0.01017929, 0.00657398, 0.01017597, 0.00657434, 0.01018185, 0.01017295, 0.01017918, 0.01018114, 0.01017743, 0.00657273, 0.01017534, 0.01017491, 0.00657773, 0.01018109],
				[0.01017665, 0.0101768, 0.00657714, 0.0101759, 0.00656631, 0.01017259, 0.01017779, 0.01017938, 0.00657658, 0.01018236, 0.00656997, 0.01018064, 0.01017251, 0.01017437, 0.0065706],
				[0.01017233, 0.0101795, 0.01017447, 0.01017058, 0.00697279, 0.00697084, 0.01017066, 0.00698027, 0.0069664, 0.0101718, 0.00696832, 0.00697924, 0.01017795, 0.01017431, 0.01017373],
				[0.01016847, 0.01016698, 0.00656584, 0.01017134, 0.00657776, 0.01017129, 0.01017305, 0.00657981, 0.01017648, 0.01017069, 0.00656734, 0.01016605, 0.01016577, 0.01016411, 0.00657812]]
				#print 't:',self.time,'neuronIndex',neuronIndex,'main calc resistanceCalc()', testR

				self.testR = testR

				#testId[neuronIndex], self.testIdSpikeFired[neuronIndex], self.testIdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent2(neuronIndex, testId, self.testIdSpikeFired, self.testIdRefractoryPeriod, testR, testW, self.testSpiketimes, True)
				
				for indexOfDend in range(dictionaryLongitude):
					# TODO see if some of these do not need to be recomputed every time
					testDend[indexOfDend].w = testW[indexOfDend]*volt
					testDend[indexOfDend].tau = testTauD[indexOfDend]*ms
					testDend[indexOfDend].r = testR[indexOfDend]*volt	

			def OutputEvaluationResults(self):
				print 'Final Weights\n0',dend[0].w[:]
				print '1',dend[1].w[:]
				print '2',dend[2].w[:]
				print '3',dend[3].w[:]
				print 'Final Res\n0',dend[0].r[:]
				print '1',dend[1].r[:]
				print '2',dend[2].r[:]
				print '3',dend[3].r[:]		
				print '\n'
				print '+++ Results +++'
				print 'Spike results: TP:\t',self.truePositiveSpikeResults,'\tFP:\t',self.falsePositiveSpikeResults,'\tTN:\t',self.trueNegativeSpikeResults,'\tFN:\t',self.falseNegativeSpikeResults
				print 'totalSpikeIntervalsTested:\t',self.totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
				print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(self.truePositiveSpikeResults, '.1f'))/Decimal(format(self.totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
				print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',(Decimal(format(self.truePositiveSpikeResults, '.1f'))+Decimal(format(self.trueNegativeSpikeResults, '.1f')))/(Decimal(format(self.totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
				print '+++++++++++++++'

			intitializeParameters()

		#S = Synapses(testADDS, testADDS, pre='vTest -= (vTest*.9)')
		#S.connect('i != j') # all-to-all but no self-connections

		ADDS = ADDSNeuronModel(4)
		testADDS = testADDSNeuronModel(4)	

		def lateralInhibition(neuronIndex):
			for inhibitionNeuronIndex in range(len(testADDS.vTest)):
				if inhibitionNeuronIndex != neuronIndex:
					testADDS.vTest[inhibitionNeuronIndex] -= testADDS.vTest[neuronIndex]*0.5
					#a = 1

		M = SpikeMonitor(ADDS)
		self.M = M # for ipython compatibility
		testM = SpikeMonitor(testADDS)
		self.testM = testM
		Mv = StateMonitor(ADDS, 'V', record=True)
		testMv = StateMonitor(testADDS, 'v', record=True)
		MDendI = StateMonitor(ADDS, 'DendI', record=True)
		MSynI = StateMonitor(ADDS, 'SynI', record=True)
		UmM3 = StateMonitor(ADDS, 'v2', record=True)
		self.UmM3 = UmM3 # for ipython compatibility
		testUmM3 = StateMonitor(testADDS, 'v2', record=True)
		self.testUmM3 = testUmM3
		dendSpikeM = [None]*dictionaryLongitude
		dendVoltageM = [None]*dictionaryLongitude
		testDendSpikeM = [None]*dictionaryLongitude
		testDendVoltageM = [None]*dictionaryLongitude
		for firstLayerIndex in range(dictionaryLongitude):
			dendSpikeM[firstLayerIndex] = SpikeMonitor(dend[firstLayerIndex])
			dendVoltageM[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'v', record=True)
			testDendSpikeM[firstLayerIndex] = SpikeMonitor(dend[firstLayerIndex])
			testDendVoltageM[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'v', record=True)

			neuralnet.add(dend[firstLayerIndex])
			neuralnet.add(dendSpikeM[firstLayerIndex])
			neuralnet.add(dendVoltageM[firstLayerIndex])
			neuralnet.add(testDend[firstLayerIndex])
			neuralnet.add(testDendSpikeM[firstLayerIndex])
			neuralnet.add(testDendVoltageM[firstLayerIndex])			

		neuralnet.add(ADDS)
		neuralnet.add(testADDS)
		neuralnet.add(M)
		neuralnet.add(testM)
		neuralnet.add(testMv)
		neuralnet.add(UmM3)
		neuralnet.add(testUmM3)
		neuralnet.run(129*ms,report='text')

		testADDS.OutputEvaluationResults()

		neuronToPlot = 1
		subplot(221)
		plot(M.t/ms, M.i, '.')
		#plot(M.t/ms, M.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(222)
		plot(testM.t/ms, testM.i, '.')
		#plot(testM.t/ms, testM.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(223)
		#plot(UmM3.t, UmM3.v2.T/mV)
		plot(UmM3.t, UmM3.v2.T/mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		subplot(224)
		plot(testUmM3.t, testUmM3.v2.T/mV)		
		#plot(testUmM3.t, testUmM3.v2Test.T/mV)
		#plot(testMv.t, testMv.vTest.T/mV)
		#plot(dendVoltageM[0].t, dendVoltageM[0].v.T/mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		show()

	def __init__(self):
		self.run_model()

def main():
	run_gupta_paper = gupta_paper()

if  __name__ =='__main__':main()
