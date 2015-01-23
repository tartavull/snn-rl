from architecture_further_formulas import *

class gupta_paper:
	neuralnet = Network()
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	testSpiketimes = spiketimes
	print 'st',spiketimes[0:100]
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	neuronIndex = 0
	time = 0.0000
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
	currentEpoch = 0
	refractoryPointCounter = 0.0
	refractoryPointSwitch = False
	neuronIndexCounter = 0.0
	neuronIndexSwitch = False
	timing = []
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
	testUmSpikeFired = np.array([False]*dictionaryLongitude); 
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

	print 'initial Weights\n',W

	def run_model(self):
		neuralnet = self.neuralnet
		dictionary = self.dictionary
		spiketimes = self.spiketimes
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
			testUmSpikeFired : boolean	
			lateralInhibActive : boolean	
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
		testADDS = NeuronGroup(N=4, model=testEqs,threshold='vTest>.01*mV', reset='vTest=-0.002 * mV; dvTest=0; IdVZTest = 0*mV;v2Test=10*mV;testUmSpikeFired=True;lateralInhibActive=True',refractory=10*ms,clock=testMyclock)'''

		class ADDSNeuronModel(NeuronGroup, gupta_paper): 
			neuronIndex = self.neuronIndex
			refractoryPointCounter = self.refractoryPointCounter
			def __init__(self, params): 
				myclock=Clock(dt=1*ms)
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV',refractory=10*ms,clock=myclock)
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.002*mV', reset='v=-0.002 * mV',refractory=0.3*ms) 
				@network_operation 
				def additionToNetwork(): 
					neuronIndex = self.neuronIndex
					spikeIntervalCounter = (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10
					print 'neuronIndex', self.neuronIndex, 'self.time', self.time,'ADDS.t', ADDS.t, 'ADDS.v', ADDS.v[neuronIndex]
					#print 'neuronIndex', neuronIndex, 'self.time', self.time, 'refractoryPointCounter',self.refractoryPointCounter,'ADDS.t', ADDS.t, 'testADDS.vTest', testADDS.vTest
					print 'neuronIndex', self.neuronIndex, 'self.time', self.time, 'ADDS.t', ADDS.t, 'ADDS.IdVZ', ADDS.IdVZ[neuronIndex]
					#print 'neuronIndex', neuronIndex, 'self.time', self.time, 'ADDS.t', ADDS.t, 'testADDS.IdVTest', testADDS.IdVTest

					#returnUm(self)

					def timePeriodAndRefractoryCalcs():
						neuronIndex = self.neuronIndex
						refractoryPointCounter = self.refractoryPointCounter
						refractoryPointSwitch = self.refractoryPointSwitch
						neuronIndexCounter = self.neuronIndexCounter
						neuronIndexSwitch = self.neuronIndexSwitch
						time = self.time
						
						# changed time step time for testing
						#timeStepInterval = 0.01
						#timeStepInterval = 0.002
						timeStepInterval = 0.001
						time = time + timeStepInterval
						self.timing.append(time)
						refractoryPointCounter = refractoryPointCounter + timeStepInterval
						# if statement below corrects for offset of spiketimes starting at .1 sec.
						if time >= .1:
							neuronIndexCounter = neuronIndexCounter + timeStepInterval

						# At the end of each spike time interval refractory period is turned off and weight changes
						# are calculated.  Refractory turning off here * is not correct * because it can occur for
						# less than the set refractory period.  I am just using it for the time being for testing.
						refractoryPointSwitch = False
						if refractoryPointCounter >= spikeIntervalUnformatted:
							#refractoryPointCounter = 0.009
							refractoryPointCounter = 0.000
							refractoryPointSwitch = True
							self.lastSpikeInterval = self.time
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
						self.refractoryPointSwitch = refractoryPointSwitch
						self.neuronIndexCounter = neuronIndexCounter			
						self.neuronIndexSwitch = neuronIndexSwitch

					def diracCalc(IDend, evaluationActive):
						# This method only calculates dirac
						time = self.time
						e = math.e
						#print ':::range(len(IDend[neuronIndex][:])\t','self.time\t',self.time,'\tneuronIndex\t',neuronIndex,'\t',IDend,'\t',IDend[neuronIndex]

						# normalize t to count just time within spike interval.  Formula from article just uses t in a way
						# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
						tNorm = time - (floor(time/spikeIntervalUnformatted) * spikeIntervalUnformatted) 

						dendGroup = [None]*len(dend[neuronIndex][:])
						for IdIndex in range(len(dend[neuronIndex][:])):
							tauDen = tauD[neuronIndex][IdIndex]
							#r = R[neuronIndex][IdIndex]
							#w = W[neuronIndex][IdIndex]

							# set tpresyn initially too far to trigger dirac
							#tPreSyn = -t - 1
							tPreSyn = time + 1
							for presynInd in range(shape(spiketimes)[0]):
								comparedSpikeTime = spiketimes[presynInd][1]

								# this is true? -> commented out logic below due to spiketimes perhaps not being in sorted order, a relevant entry could occur after an irrelevant one
								#if (comparedSpikeTime-spikeIntervalUnformatted) > t:
								# spike times are in order of earliest to latest times grouped with input pixel number.  therefore the below cutoff in the loop is fine.
								#if comparedSpikeTime > self.lastSpikeInterval:
								if comparedSpikeTime > (self.lastSpikeInterval + spikeIntervalUnformatted):
									break

								#print ':::',comparedSpikeTime,self.lastSpikeInterval,spikeIntervalUnformatted
								##print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'spiketimes[presynInd][0] == IdIndex',(spiketimes[presynInd][0] == IdIndex),'IdIndex',IdIndex,'comparedSpikeTime > self.lastSpikeInterval',(comparedSpikeTime > self.lastSpikeInterval),'comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted)',(comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted))
								'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
								# 0.0 to 300 ms, as the input is intended
								# changed <= to < below to see if that works better
								if spiketimes[presynInd][0] == IdIndex and (comparedSpikeTime-spikeIntervalUnformatted) < t:'''
								# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
								#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted):
								#print 'self.time',self.time,'neuronIndex',neuronIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0],'comparedSpikeTime',comparedSpikeTime#,'self.lastSpikeInterval',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'IdIndex',IdIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0]
								#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - (spikeIntervalUnformatted*2)) and comparedSpikeTime <= (self.lastSpikeInterval - spikeIntervalUnformatted):
								#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted) and comparedSpikeTime <= self.lastSpikeInterval:
								if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > self.lastSpikeInterval and comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted):
									## adding below if statement to avoid constant input stimulus in evaluation, that seem more accurate to how the article did it's evaulation.
									#print 'tNorm',tNorm
									tNorm2 = time - (floor((time/.001)*.01) * .1)
									print 'tNorm2',tNorm2
									if tNorm2 < .02 or evaluationActive == False:
										tPreSyn = comparedSpikeTime
										print '++++++++++'

							#Id2 = IDend[neuronIndex].v[IdIndex]
							#Dt = time - tPreSyn
							Dt = Decimal(format(time, '.8f')) - Decimal(format(tPreSyn, '.8f'))
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
							tauDen = tauDen * .001
							
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
								dendGroup[IdIndex] = float(DiracFun)*volt						

							#print 'neuronIndex',neuronIndex,'IdIndex',IdIndex,'dirac\t',DiracFun
						return dendGroup

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
						dend[neuronIndex].dirac = diracCalc(dend, False)
						print 'dend.dirac',dend[0].dirac,dend[1].dirac,dend[2].dirac,dend[3].dirac
						print 'dend.v',dend[0].v,dend[1].v,dend[2].v,dend[3].v
						print 'sum(dend[neuronIndex].v[:])',sum(dend[neuronIndex].v[:])
						'''for i in range(4):
							if i != neuronIndex:
								dend[i].dirac = [0*volt]*15'''

						# Calculate tauD
						tauD = self.tauD
						tauD = tauDCalc(neuronIndex, dend, W)
						#self.tauD = tauD # this line is just for resutls reporting and can be removed later
						print 'tauD',tauD

						# Calculate resistance
						R = self.R
						R = resistanceCalc(neuronIndex, dend, R)
						#self.R = R # this line is just for resutls reporting and can be removed later
						print 'R',R

						# Record sum of dend voltage
						#ADDS.IdVZ[neuronIndex] = sum(dend[neuronIndex].v[:])#0.0375*mV 
						for i in range(4):
							ADDS.IdVZ[i] = sum(dend[i].v[:])
						'''for i in range(4):
							if i != neuronIndex:
								ADDS.IdVZ[i] = 0*volt'''
						print '::ADDS.IdVZ::\t', ADDS.IdVZ
						print 'ADDS.v2',ADDS.v2

						timePeriodAndRefractoryCalcs()

					mainSimulationCalcs()

				self.contained_objects.append(additionToNetwork)		

		class testADDSNeuronModel(NeuronGroup, gupta_paper): 
			def __init__(self, params): 
				testMyclock=Clock(dt=1*ms)
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV;testUmSpikeFired=True;lateralInhibActive=True',refractory=10*ms,clock=testMyclock)
				#NeuronGroup.__init__(self, N=4, model=testEqs,threshold='vTest>.002*mV', reset='vTest=-0.002 * mV',refractory=0.3*ms) 
				#NeuronGroup.__init__(self, N=15, model=dendriteEqs,threshold='v>.002*mV', reset='v=-0.002 * mV',refractory=0.3*ms) 
				@network_operation 
				def additionToNetwork(): 
					spikeIntervalCounter = (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10

					for vCheck in range(len(testADDS.v2Test)):
						if testADDS.v2[vCheck] == 10*mV:
							testADDS.v2[vCheck] = 0*mV					

					# classifier performance test
					#evaluateClassifier()

					# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
					# the +1 in (self.timeotalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
					if spikeIntervalCounter < (self.totalSpikeIntervals+1):
						for neuronIndex in range(dictionaryLongitude):				
							# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
							# (spikeIntervalCounter-1) is to correct for refractoryPointSwitch occuring after spikeInterval it addresses.
							if self.refractoryPointSwitch == true and (spikeIntervalCounter > 0):
								if self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter-1] == False:
									if (self.correctSpikes[neuronIndex][(spikeIntervalCounter-1)] == 1):
										self.falseNegativeSpikeResults = self.falseNegativeSpikeResults + 1		
									else:
										self.trueNegativeSpikeResults = self.trueNegativeSpikeResults + 1	

								if (testADDS.testUmSpikeFired[neuronIndex] == True) and (spikeIntervalCounter < self.totalSpikeIntervals):
									if (self.correctSpikes[neuronIndex][spikeIntervalCounter] == 1):
										self.truePositiveSpikeResults = self.truePositiveSpikeResults + 1	
									else:
										self.falsePositiveSpikeResults = self.falsePositiveSpikeResults + 1	
									self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = True	
									testADDS.testUmSpikeFired[neuronIndex] = False
						self.contained_objects.append(additionToNetwork)


					def evaluateClassifier():
						#tNorm = self.time - (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
						tNorm = time - (floor((time/.001)*.01) * .1)

						'''testW = [[1, 0, 1,0, 1, 0,0, 0, 0,0, 1, 0,0, 1, 0],
						[0, 0, 1,0, 1, 0,0, 0, 1,0, 1, 0,0, 0, 1],
						[0, 0, 0,0, 1, 1,0, 1, 1,0, 1, 1,0, 0, 0],
						[0, 0, 1,0, 1, 0,0, 1, 0,0, 1, 0,0, 0, 1]]'''
						testW = [[0.70080412, 0.00324985, 0.70223351, 0.00421194, 0.70073588, 0.00253654, 0.00241304, 0.00386221, 0.00446284, 0.0045989, 0.70021551, 0.00332156, 0.00453394, 0.70026736, 0.00364166],
						[0.00312473, 0.00303151, 0.70134693, 0.00250971, 0.70082969, 0.00264048, 0.00426394, 0.00245731, 0.70103305, 0.00272881, 0.70223473, 0.00287218, 0.00364664, 0.00237047, 0.70049339],
						[0.00337244, 0.00489955, 0.00332976, 0.00476442, 0.66939278, 0.66919182, 0.00383946, 0.66793638, 0.6669667, 0.00301663, 0.6670524, 0.66844534, 0.00401786, 0.00284326, 0.00354638],
						[0.00381406, 0.00498837, 0.701597, 0.00521621, 0.70256855, 0.00533906, 0.00412508, 0.7023335, 0.00631239, 0.0043698, 0.70042287, 0.00589082, 0.00326269, 0.00518741, 0.70275034]]

						# for each charactor input below test the classfication results
						for neuronIndex in range(dictionaryLongitude):
							# Calculate tauD
							testTauD = self.testTauD
							testTauD = tauDCalc(neuronIndex, testTauD, testW)
							#print 't:',self.time,'neuronIndex',neuronIndex,'main calc tauDCalc()', testTauD

							testTauD = [10.37748461, 29.90900423, 10.33746184, 29.88206559, 10.37939543, 29.928977, 29.932435, 29.89185801, 29.87504059, 29.87123073, 10.39396579, 29.90699629, 29.87304979, 10.39251402, 29.89803348]
							[29.91250761, 29.91511768, 10.36228604, 29.92972817, 10.37676875, 29.92606665, 29.8806098, 29.93119539, 10.37107451, 29.92359325, 10.33742747, 29.91957909, 29.89789419, 29.93362677, 10.38618516]
							[29.89507952, 29.84756941, 29.89640735, 29.85177368, 12.28555811, 12.29181016, 29.88055016, 12.3308681, 12.3610359, 29.90614938, 12.35836972, 12.31503375, 29.87499997, 29.911543, 29.88966829]
							[29.89320623, 29.86032567, 10.35528394, 29.85394614, 10.32808049, 29.85050638, 29.88449768, 10.33466202, 29.82325308, 29.87764556, 10.38815977, 29.83505713, 29.90864456, 29.85475254, 10.32299042]

							# Calculate resistance
							testR = self.testR
							testR = resistanceCalc(neuronIndex, testTauD, testR)

							#Sample output for testing
							'''testR = [[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
							[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
							[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
							[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463]]'''
							testR = [[0.70080412, 0.00324985, 0.70223351, 0.00421194, 0.70073588, 0.00253654, 0.00241304, 0.00386221, 0.00446284, 0.0045989, 0.70021551, 0.00332156, 0.00453394, 0.70026736, 0.00364166],
							[0.00312473, 0.00303151, 0.70134693, 0.00250971, 0.70082969, 0.00264048, 0.00426394, 0.00245731, 0.70103305, 0.00272881, 0.70223473, 0.00287218, 0.00364664, 0.00237047, 0.70049339],
							[0.00337244, 0.00489955, 0.00332976, 0.00476442, 0.66939278, 0.66919182, 0.00383946, 0.66793638, 0.6669667, 0.00301663, 0.6670524, 0.66844534, 0.00401786, 0.00284326, 0.00354638],
							[0.00381406, 0.00498837, 0.701597, 0.00521621, 0.70256855, 0.00533906, 0.00412508, 0.7023335, 0.00631239, 0.0043698, 0.70042287, 0.00589082, 0.00326269, 0.00518741, 0.70275034]]
							#print 't:',self.time,'neuronIndex',neuronIndex,'main calc resistanceCalc()', testR

							#testId[neuronIndex], self.testIdSpikeFired[neuronIndex], self.testIdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent2(neuronIndex, testId, self.testIdSpikeFired, self.testIdRefractoryPeriod, testR, testW, self.testSpiketimes, True)
							
							testDend[neuronIndex].dirac = ADDS.network_operation.additionToNetwork.diracCalc(testDend, False)

							for indexOfDend in range(dictionaryLongitude):
								# TODO see if some of these do not need to be recomputed every time
								testDend[indexOfDend].w = testW[indexOfDend]
								testDend[indexOfDend].tau = testTauD[indexOfDend]
								testDend[indexOfDend].r = testR[indexOfDend]
								testADDS.IdVZ[indexOfDend] = sum(testDend[indexOfDend].v[:])

							#print 'test neuronIndex2',neuronIndex,'self.time',self.time,'testADDS.t',testADDS.t,'testId[neuronIndex]',testId[neuronIndex],'testADDS.v2Test',testADDS.v2Test
							#testIs = None
							#testADDS.IdVTest[neuronIndex] = totalSomaMembranePotential(neuronIndex, testId, testIs, tNorm)		

		#S = Synapses(testADDS, testADDS, pre='vTest -= (vTest*.9)')
		#S.connect('i != j') # all-to-all but no self-connections

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

		ADDS = ADDSNeuronModel(4)
		testADDS = testADDSNeuronModel(4)
		dend = [None]*dictionaryLongitude
		testDend = [None]*dictionaryLongitude
		for firstLayerIndex in range(dictionaryLongitude):
			dend[firstLayerIndex] = DendriteNeuronModel(15)
			testDend[firstLayerIndex] = DendriteNeuronModel(15)

		def returnUm(self):
			# Main logic of the neuron simulation model is performed here
			# Modules are performed below resulting in voltage for the neurons being calculated

			dictionary = self.dictionary
			neuronIndex = self.neuronIndex
			time = self.time
			tNorm = time - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
			spiketimes = self.spiketimes
			
			# Calculate tauD
			tauD = self.tauD
			tauD = tauDCalc(neuronIndex, tauD, W)
			#print 'main calc tauDCalc()', tauD

			# Calculate resistance
			R = self.R
			R = resistanceCalc(neuronIndex, tauD, R)
			#print 'main calc resistanceCalc()', R

			'''#Dedritic total post-synaptic current
			Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, Id, self.IdSpikeFired, self.IdRefractoryPeriod, R, W, spiketimes, False)
			#print 'neuronIndex',neuronIndex,'self.time',self.time,'ADDS.t',ADDS.t,'Id[neuronIndex]',Id[neuronIndex]
			# Direct to soma 

			### Soma membrane potential ###
			ADDS.IdVZ[neuronIndex] = totalSomaMembranePotential(neuronIndex, Id, Is, tNorm)

			# Refractory and spike evaluation: Is
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])			'''

			timePeriodAndRefractoryCalcs()

			return Um[neuronIndex] * mV;	

		def refractoryPeriodEvaluation(Voltage, SpikeFired, RefractoryPeriod):
			# Check for spike
			# when identifying a spike, refractoryPeriodEndPoint is used as an alternative condition to refractoryPeriod being
			# over to allow a spike right at that end point

			# Refractory period
			#refractoryPeriodEndPoint = (self.time - 0.001) % spikeIntervalUnformatted <= 0.00001			
			if SpikeFired == True and self.refractoryPointSwitch == False:
				RefractoryPeriod = True
				#print 'refrac activated!'
			elif RefractoryPeriod == True and self.refractoryPointSwitch == True:
				# 0.001 added above for time step compatibility.   0.00001 instead of 0.0 used for some kind of 
				# modulo computational result offset in python which was producing a value from the modulo calc
				# just slightly over 0
				RefractoryPeriod = False
				#print 'refrac over!'
			#print 'self.refractoryPeriod == True and refractoryPeriodEndPoint', self.refractoryPeriod == True, refractoryPeriodEndPoint, (self.refractoryPeriod == True and refractoryPeriodEndPoint), self.time % (spikeIntervalUnformatted + 0.001) <= 0.00001
			if RefractoryPeriod == True:
				Voltage = Ureset			

			return [Voltage, SpikeFired, RefractoryPeriod]

		def spikeFiredEvaluation(Voltage, SpikeFired, RefractoryPeriod):			
			# Spike firing detection positioned here to allow Ureset value to occur on next timestep.
			# This allows spike to be recorded in plot instead of value immediately going to Ureset.
			SpikeFired = False
			if Voltage >= ActionPotentialThreshold and (RefractoryPeriod == False or refractoryPeriodEndPoint):
				#Um[neuronIndex] = Ureset
				#print 'Voltage',Voltage
				SpikeFired = True	

			return [Voltage, SpikeFired, RefractoryPeriod]

		def dentritePostSynapticCurrent(neuronIndex, IDend, IDendSpikes, IDendRefract, R, W, spiketimes, evaluationActive):
			# Solving for Idend in the formula in the article yeilded the below equation
			time = self.time
			e = math.e
			#print ':::range(len(IDend[neuronIndex][:])\t','self.time\t',self.time,'\tneuronIndex\t',neuronIndex,'\t',IDend,'\t',IDend[neuronIndex]

			# normalize t to count just time within spike interval.  Formula from article just uses t in a way
			# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
			tNorm = time - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 

			for IdIndex in range(len(IDend[neuronIndex][:])):
				tauDen = tauD[neuronIndex][IdIndex]
				r = R[neuronIndex][IdIndex]
				w = W[neuronIndex][IdIndex]

				# set tpresyn initially too far to trigger dirac
				#tPreSyn = -t - 1
				tPreSyn = time + 1
				for presynInd in range(shape(spiketimes)[0]):
					comparedSpikeTime = spiketimes[presynInd][1]
					#print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval
					'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
					# 0.0 to 300 ms, as the input is intended
					# changed <= to < below to see if that works better
					if spiketimes[presynInd][0] == IdIndex and (comparedSpikeTime-spikeIntervalUnformatted) < t:'''
					# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
					#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted):
					#print 'self.time',self.time,'neuronIndex',neuronIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0],'comparedSpikeTime',comparedSpikeTime#,'self.lastSpikeInterval',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'IdIndex',IdIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0]
					if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime >= (self.lastSpikeInterval - spikeIntervalUnformatted) and comparedSpikeTime <= self.lastSpikeInterval:
						## adding below if statement to avoid constant input stimulus in evaluation, that seem more accurate to how the article did it's evaulation.
						#print 'tNorm',tNorm
						tNorm2 = time - (floor((t/.001)*.01) * .1)
						print 'tNorm2',tNorm2
						if tNorm2 < .02 or evaluationActive == False:
							tPreSyn = comparedSpikeTime

					# this is true? -> commented out logic below due to spiketimes perhaps not being in sorted order, a relevant entry could occur after an irrelevant one
					#if (comparedSpikeTime-spikeIntervalUnformatted) > t:
					# spike times are in order of earliest to latest times grouped with input pixel number.  therefore the below cutoff in the loop is fine.
					if comparedSpikeTime > self.lastSpikeInterval:
						break

				#print 'tPreSyn000:',tPreSyn
				Id2 = IDend[neuronIndex][IdIndex]
				Dt = time - tPreSyn
				# self.time == lines below are a workaround for initialization values
				if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
					Dt = 1.0
				# dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
				# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
				# though
				# simplify dirac for testing
				#DiracFun = 0
				DiracFun = 1/Dt

				#SpikeModCoeff = (r*w*DiracFun)

				# dirac test
				# t in dirac forumula means curent time or last spike time?		
				#if (t > -(Dt/2) and t < (Dt/2)):
				#if Dt <= spikeIntervalUnformatted:
				SpikeModCoeff = (r*w)
				if evaluationActive == True: SpikeModCoeff *= 0.018#0.018#0.015#0.003#0.01#0.0206325
				
				#print '***ni:\t',neuronIndex,'\tIdIndex\t',IdIndex,'\t***w:\t',w

				# correct for scaling
				tauDen = tauDen * .001
				
				#if -Dt/2<(tNorm or t?)<Dt/2:
				#if Dt <= 0.0:
				if Dt >= 0.0:
					if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
						DiracFun = 1.00
					#SpikeModCoeff = .011
					#DiracFun = 0.66
					SpikeModCoeff = (SpikeModCoeff*DiracFun)
					#tauDen = .03	
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
					if evaluationActive == True: print 'test IDend2:\tneuronIndex',neuronIndex,'IdIndex',IdIndex,'self.time',self.time,'IDend[neuronIndex]',IDend[neuronIndex],'SpikeModCoeff',SpikeModCoeff
				else:
					DiracFun = 0
					SpikeModCoeff = (SpikeModCoeff*DiracFun)															
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
					#print 'IDend2:\t',self.time,neuronIndex,IDend[neuronIndex]
				
				#print 'part1: \t',IDend
				#print 'factors\tSpikeModCoeff\t',SpikeModCoeff,'self.time',self.time,'neuronIndex',neuronIndex,'\tId2\t',Id2,'\ttPreSyn\t',tPreSyn,'\tDiracFun\t',DiracFun,'Dt <= 0.0',(Dt <= 0.0),'Dt',Dt,'\tneuronIndex\t',neuronIndex,'\tIdIndex\t',IdIndex,'-(SpikeModCoeff - Id2)',-(SpikeModCoeff - Id2),'-(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen))',-(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen))

				# Refractory and spike evaluation: Id
				#print 'part2: \t',refractResults
				#spikeEvalResults = spikeFiredEvaluation(refractResults[0], refractResults[1], refractResults[2])
				#print 'part3: \t',spikeEvalResults
				spikeEvalResults = IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex]
				if self.time > (self.lastSpikeInterval+0.02):
					refractResults = refractoryPeriodEvaluation(IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex])
					spikeEvalResults = spikeFiredEvaluation(refractResults[0], refractResults[1], refractResults[2])
				
				IDend[neuronIndex][IdIndex] = spikeEvalResults[0]
				#print 'IDend3:\t',self.time,neuronIndex,IDend[neuronIndex]
				IDendSpikes[neuronIndex][IdIndex] = spikeEvalResults[1]
				IDendRefract[neuronIndex][IdIndex] = spikeEvalResults[2]

				#print 'part1_2: \t', IDend
				
			return [IDend[neuronIndex], IDendSpikes[neuronIndex], IDendRefract[neuronIndex]]

		def dentritePostSynapticCurrent2(neuronIndex, modelOfNeurons, R, W, spiketimes, evaluationActive):
			modelOfNeurons.r = 80*mV # change later
			modelOfNeurons.w = 1
			modelOfNeurons.dirac = 1
			modelOfNeurons.tau = 30*ms # change later

		def somaPostSynapticCurrent(neuronIndex, ISoma, ISomaSpikes, ISomaRefract, R, W, spiketimes):
			'''### Synapse directly to soma ###
			# Solving for Id in the formula in the article yeilded the below equation
			# To calculate the DiracWeightedSum the spike times with the dirac function applied are multipled by the synapse weight 
			# and summed then divided by the number of synapses for the neuron 
			DiracWeightedSum = 0
			for DiracIndex in range(len(Is)):
				#tPreSyn = spiketimes[DiracIndex + (neuronIndex * len(Is))][0]
				tPreSyn = spiketimes[DiracIndex + (neuronIndex * len(Is)) + (self.currentEpoch*(dictionaryLongitude*numberOfPixels))][0]
				Dt = t - tPreSyn
				DiracFun = 1/Dt
				# dirac test  # TODO: not sure dirac is implemented correctly here
				if (t > -(Dt/2) and t < (Dt/2)):
					DiracFunctionWithSpikeTimes = DiracFun # (-t+tPreSyn)/tauDen
				else:
					DiracFunctionWithSpikeTimes = 0
				DiracWeightedSum = DiracWeightedSum + W[neuronIndex][WIndex] * DiracFunctionWithSpikeTimes
			DiracWeightedSum = DiracWeightedSum / len(Id)

			Is2 = Is[neuronIndex]

			#Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-t/tauS)) + DiracWeightedSum
			Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-tNorm/tauS)) + DiracWeightedSum'''

			return ['stub', 'stub', 'stub']

		def totalSomaMembranePotential(neuronIndex, IDend, ISoma, tNorm):
			### Soma membrane potential.  Now this method is just used for dendrite synapse output aggregation ###
			SummedDendriteGroup = sum(IDend[neuronIndex])*volt#*10#*10 added 12/29/14 for scaling adjustment test # * 1 # the * 1000 is for a scaling adjustment
			# below removed until later needed
			#SynapseToSoma = ISoma[neuronIndex]
			#print 'self.epochIndex\t',self.epochIndex,'\tneuronIndex\t',neuronIndex,'\tSummedDendriteGroup\t',SummedDendriteGroup,'\tIDend\t',IDend
			return SummedDendriteGroup

		'''def WeightChangeCalculation():
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
			##			'''

		def evaluateClassifier():
			tNorm = self.time - (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted) 

			testR = [[1] * numberOfPixels]*dictionaryLongitude
			'''testW = [[1, 0, 1,0, 1, 0,0, 0, 0,0, 1, 0,0, 1, 0],
			[0, 0, 1,0, 1, 0,0, 0, 1,0, 1, 0,0, 0, 1],
			[0, 0, 0,0, 1, 1,0, 1, 1,0, 1, 1,0, 0, 0],
			[0, 0, 1,0, 1, 0,0, 1, 0,0, 1, 0,0, 0, 1]]'''
			'''testW = [[ 0.69607004, 0.00413144, 0.91697601, 0.00360609, 0.91761815, 0.00445725, 0.00340361, 0.22605207, 0.0039314, 0.00451545, 0.91885474, 0.00333926, 0.00245632, 0.69608306, 0.22541202],
			[ 0.29794995, 0.0025975, 0.99158123, 0.00318387, 0.99296564, 0.00275211, 0.00305578, 0.0027208, 0.69734084, 0.00306026, 0.99202054, 0.00387052, 0.00440005, 0.29888318, 0.69716381],
			[ 0.0048261, 0.00419521, 0.33389914, 0.00392254, 0.96962399, 0.64231264, 0.00371715, 0.64278451, 0.97156676, 0.00332437, 0.96967301, 0.64113085, 0.00312674, 0.00290903, 0.33244029],
			[ 0.00562191, 0.00532413, 0.66298468, 0.00350251, 0.9567765, 0.29948389, 0.00621025, 0.9585158, 0.30201652, 0.00402158, 0.95678115, 0.30231021, 0.00520207, 0.00577489, 0.66244477]]'''
			testW = [[0.69882415, 0.00375697, 0.89738531, 0.00382883, 0.89708958, 0.00252481, 0.00245958, 0.20249973, 0.00421835, 0.00289073, 0.89693852, 0.0036611, 0.0044831, 0.69794798, 0.20160954],
			[0.20602601, 0.003453, 0.90187063, 0.00367539, 0.9013504, 0.00233126, 0.00431821, 0.00276384, 0.69783304, 0.00371228, 0.90245003, 0.00262979, 0.00285232, 0.20685335, 0.69778261],
			[0.00482211, 0.00497251, 0.22888016, 0.00452582, 0.86717426, 0.64120109, 0.00297225, 0.64273965, 0.86776953, 0.00508325, 0.86628416, 0.64287156, 0.00407194, 0.00322563, 0.2305038],
			[0.00341097, 0.00332899, 0.65583, 0.00590996, 0.85431181, 0.20267192, 0.00589263, 0.85207654, 0.20215162, 0.00352097, 0.85206755, 0.19999237, 0.00586857, 0.00587151, 0.6560445]]

			# for each charactor input below test the classfication results
			for neuronIndex in range(dictionaryLongitude):
				# Calculate tauD
				testTauD = self.testTauD
				testTauD = tauDCalc(neuronIndex, testTauD, testW)
				#print 't:',self.time,'neuronIndex',neuronIndex,'main calc tauDCalc()', testTauD

				# Calculate resistance
				testR = self.testR
				testR = resistanceCalc(neuronIndex, testTauD, testR)

				#Sample output for testing
				'''testR = [[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463],
				[0.01016803758319553, 0.010171374818041036, 0.0083665317732988555, 0.01017281704840039, 0.0048279096384584727, 0.0071335782968518762, 0.010173903373815114, 0.0071306637775638056, 0.0048108924488047102, 0.010175980702300289, 0.0048274813185357644, 0.0071408747375879501, 0.010177025843651875, 0.010178177167130422, 0.0083748475057425463]]'''
				testR = [[0.010168058717316278, 0.010167263109752002, 0.0089565412798334777, 0.010169626034713717, 0.0056387671561491737, 0.0071404412116810954, 0.010177842862684333, 0.0071309409182563821, 0.0056344099035213341, 0.010166677236855375, 0.0056452767455983871, 0.0071301260300992786, 0.010172026813988511, 0.010176502893209517, 0.0089475431290829757],
				[0.010168058717316278, 0.010167263109752002, 0.0089565412798334777, 0.010169626034713717, 0.0056387671561491737, 0.0071404412116810954, 0.010177842862684333, 0.0071309409182563821, 0.0056344099035213341, 0.010166677236855375, 0.0056452767455983871, 0.0071301260300992786, 0.010172026813988511, 0.010176502893209517, 0.0089475431290829757],
				[0.010168058717316278, 0.010167263109752002, 0.0089565412798334777, 0.010169626034713717, 0.0056387671561491737, 0.0071404412116810954, 0.010177842862684333, 0.0071309409182563821, 0.0056344099035213341, 0.010166677236855375, 0.0056452767455983871, 0.0071301260300992786, 0.010172026813988511, 0.010176502893209517, 0.0089475431290829757],
				[0.010168058717316278, 0.010167263109752002, 0.0089565412798334777, 0.010169626034713717, 0.0056387671561491737, 0.0071404412116810954, 0.010177842862684333, 0.0071309409182563821, 0.0056344099035213341, 0.010166677236855375, 0.0056452767455983871, 0.0071301260300992786, 0.010172026813988511, 0.010176502893209517, 0.0089475431290829757]]

				#print 't:',self.time,'neuronIndex',neuronIndex,'main calc resistanceCalc()', testR

				testId[neuronIndex], self.testIdSpikeFired[neuronIndex], self.testIdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent2(neuronIndex, testId, self.testIdSpikeFired, self.testIdRefractoryPeriod, testR, testW, self.testSpiketimes, True)
				print 'test neuronIndex2',neuronIndex,'self.time',self.time,'testADDS.t',testADDS.t,'testId[neuronIndex]',testId[neuronIndex],'testADDS.v2Test',testADDS.v2Test
				testIs = None
				testADDS.IdVTest[neuronIndex] = totalSomaMembranePotential(neuronIndex, testId, testIs, tNorm)
				
		def tauDCalc(neuronIndex, tau, W):
			# Weight loop
			for WIndex in range(len(W[0][:])):
				if abs(W[neuronIndex][WIndex]) <= 1:
					tau[neuronIndex][WIndex] = tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)

			return tau

		def resistanceCalc(neuronIndex, tau, R):
			tauM = self.tauM
			#Resistance loop
			for RIndex in range(len(R[0][:])):				
				if (tau[neuronIndex][RIndex]*.001 == tauM):
					# avoid a division by 0 issue
					tauM = tauM - .000001

				R[neuronIndex][RIndex] = (((tau[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm) * ((tauM / (tau[neuronIndex][RIndex]*.001) ) ** (tauM / (tauM - (tau[neuronIndex][RIndex]*.001) )))

			return R

		def old_timePeriodAndRefractoryCalcs():
			neuronIndex = self.neuronIndex
			refractoryPointCounter = self.refractoryPointCounter
			refractoryPointSwitch = self.refractoryPointSwitch
			neuronIndexCounter = self.neuronIndexCounter
			neuronIndexSwitch = self.neuronIndexSwitch
			time = self.time
			
			# changed time step time for testing
			#timeStepInterval = 0.01
			#timeStepInterval = 0.002
			timeStepInterval = 0.001
			time = time + timeStepInterval
			self.timing.append(t)
			refractoryPointCounter = refractoryPointCounter + timeStepInterval
			neuronIndexCounter = neuronIndexCounter + timeStepInterval

			# At the end of each spike time interval refractory period is turned off and weight changes
			# are calculated.  Refractory turning off here * is not correct * because it can occur for
			# less than the set refractory period.  I am just using it for the time being for testing.
			refractoryPointSwitch = False
			if refractoryPointCounter >= spikeIntervalUnformatted:
				#refractoryPointCounter = 0.009
				refractoryPointCounter = 0.000
				refractoryPointSwitch = True
				self.lastSpikeInterval = self.time
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

			# revision of when new neuron index is used.  new neuron every 300ms
			if neuronIndexSwitch == True:
				neuronIndex = neuronIndex + 1

			if neuronIndex == 4: 
				neuronIndex = 0

			self.neuronIndex = neuronIndex
			self.time = time
			self.refractoryPointCounter = refractoryPointCounter
			self.refractoryPointSwitch = refractoryPointSwitch
			self.neuronIndexCounter = neuronIndexCounter			
			self.neuronIndexSwitch = neuronIndexSwitch

		def OutputEvaluationResults(W, R, dictionaryLongitude):
			print 'Final Weights\n0',W[0][:]
			print '1',W[1][:]
			print '2',W[2][:]
			print '3',W[3][:]
			R[0] = dend[0].r[:]
			R[1] = dend[1].r[:]
			R[2] = dend[2].r[:]
			R[3] = dend[3].r[:]
			print 'Final Res\n0',R[0][:]
			print '1',R[1][:]
			print '2',R[2][:]
			print '3',R[3][:]		
			print '\n'
			print '+++ Results +++'
			print 'Spike results: TP:\t',self.truePositiveSpikeResults,'\tFP:\t',self.falsePositiveSpikeResults,'\tTN:\t',self.trueNegativeSpikeResults,'\tFN:\t',self.falseNegativeSpikeResults
			print 'totalSpikeIntervalsTested:\t',self.totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
			print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(self.truePositiveSpikeResults, '.1f'))/Decimal(format(self.totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
			print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',(Decimal(format(self.truePositiveSpikeResults, '.1f'))+Decimal(format(self.trueNegativeSpikeResults, '.1f')))/(Decimal(format(self.totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
			print '+++++++++++++++'

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
		for firstLayerIndex in range(dictionaryLongitude):
			dendSpikeM[firstLayerIndex] = SpikeMonitor(dend[firstLayerIndex])
			dendVoltageM[firstLayerIndex] = StateMonitor(dend[firstLayerIndex], 'v', record=True)

			neuralnet.add(dend[firstLayerIndex])
			neuralnet.add(dendSpikeM[firstLayerIndex])
			neuralnet.add(dendVoltageM[firstLayerIndex])

		neuralnet.add(ADDS)
		neuralnet.add(testADDS)
		neuralnet.add(M)
		neuralnet.add(testM)
		neuralnet.add(testMv)
		neuralnet.add(UmM3)
		neuralnet.add(testUmM3)
		neuralnet.run(200*ms,report='text')

		#totalRunTime = 21
		#run(totalRunTime*ms,threads=2, report='text')
		#run(119*ms,report='text') # Run with enough time for scoring evaluation
		#run(2000*ms,report='text')

		# Present results and logging
		if presentResults == True:
			for epochIndex in epochsToPrint:
				reportResultsAndLogging = report_results_and_logging(dictionaryLongitude, epochsToPrint, M, Mv, epochIndex, spikeIntervalUnformatted, dictionary, epochMsDuration)
				reportResultsAndLogging.presenter()	

		# Show all membrane potentials
		# if statement with epochIndex below is to filter only epochs actaully included in a run
		if displayAllNeuronMemPotentials == True:
			for neuronIdentity in range(dictionaryLongitude):
				for epochIndex in epochsToPrint:
					if epochIndex <= totalRunTime/100:
						for epochMs in range(epochMsDuration):
							print 'neuron\t', neuronIdentity, '\tepochIndex\t', epochIndex, '\tepochMs\t', epochMs, '\tMembranePotential\t', Mv[neuronIdentity][(epochIndex*epochMsDuration)+epochMs]

		if logger == True:
			outputFile = open('snnResults.txt', 'w')	
			for epochIndex in range(epochs):
				reportResultsAndLogging = report_results_and_logging(dictionaryLongitude, epochsToPrint, M, Mv, epochIndex, spikeIntervalUnformatted, dictionary, epochMsDuration)
				reportResultsAndLogging.logger(outputFile)
			outputFile.close()

		OutputEvaluationResults(W, R, dictionaryLongitude)

		neuronToPlot = 1
		subplot(111)
		#plot(M.t/ms, M.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(221)
		plot(M.t/ms, M.i, '.')
		#plot(testM.t/ms, testM.i, '.')
		#plot(dendSpikeM[0].t/ms, dendSpikeM[0].i, '.')
		subplot(223)
		#plot(UmM3.t, UmM3.v2.T/mV)
		#plot(UmM3.t, UmM3.v2.T/mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		subplot(224)
		plot(UmM3.t, UmM3.v2.T/mV)
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
