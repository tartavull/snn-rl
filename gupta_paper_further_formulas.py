from architecture_further_formulas import *

class gupta_paper:
	neuralnet = Network()
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	trainingSpikeTimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, 1, 300)
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	absRefracTime = [0]*dictionaryLongitude	
	relRefracTime = [0]*dictionaryLongitude
	timeStepInterval = 0.001	
	time = 0.0000 - timeStepInterval
	neuronIndex = 0
	lastSpikeInterval = -0.001#0.0
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
	W = W
	R = R
	testR = testR
	tauM = tauM
	currentEpoch = 0
	refractoryPointCounter = 0.0
	testRefractoryPointCounter = 0.0
	refractoryPointSwitch = False
	testRefractoryPointSwitch = False
	neuronIndexCounter = 0.0
	neuronIndexSwitch = False
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
	IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	testIdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); 
	IsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	testIsRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	UmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	testUmRefractoryPeriod = np.array([False]*dictionaryLongitude); 
	refractoryGraphingSwitch = False

	correctSpikes = np.array([[1]*(totalSpikeIntervals+1)]*dictionaryLongitude);#[[1] * totalSpikeIntervals]*dictionaryLongitude
	correctOrderOfNeurons = [2,0,1,3]
	correctSpikes[correctOrderOfNeurons[0]][0] = 0
	correctSpikes[correctOrderOfNeurons[0]][4:13] = 0
	correctSpikes[correctOrderOfNeurons[1]][0:4] = 0
	correctSpikes[correctOrderOfNeurons[1]][7:13] = 0
	correctSpikes[correctOrderOfNeurons[2]][0:7] = 0
	correctSpikes[correctOrderOfNeurons[2]][10:13] = 0
	correctSpikes[correctOrderOfNeurons[3]][0:10] = 0
    
	M = None
	testM = None
	UmM3 = None
	testUmM3 = None
	spikeIntervalCounter = 0.0
	relRefracTimeActivated = [False]*dictionaryLongitude	
	absRefracTimeActivated = [False]*dictionaryLongitude
	# Parameters for tuning
	relRefracTimeDuration = 0.03#.1#.01		
	absRefracTimeDuration = 0.02#.02#.002	
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
	voltageForInhibition = 0.0*mV
	inhibitionReduction = 1.0
	positiveWeightReinforcement = 1.0#2.0#2.0#7.0#1.0#7.0#20.0#5.0#10.0#1.0#2.0
	negativeWeightReinforcement = 2.0#5.0#3.5#7.5#30.0#30#2.0#30#7.5#15.0#30.0#8.0#6.0#1.0#50.0#10.0#50.0#10.0#3.0#0.5
	LearningRate = 0.2#0.1#1.0#0.7
	evaluateClassifier = True
	accelerateTraining = False
	runTime = 155*ms#1000*ms#300*ms#400*ms#155*ms#42*ms#135*ms#72*ms#135*ms#42*ms#100*ms#500*ms#5000*ms#100*ms#22*ms#15*ms#135*ms#15*ms#130*ms

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
			refractoryPointCounter = self.refractoryPointCounter
			calcDirac2 = None
			generalClockDt = self.generalClockDt
			def __init__(self, params):
				NeuronGroup.__init__(self, N=4, model=eqs,threshold='v>10*mV', reset='v=-0.002 * mV; dv=0; v2=10*mV;UmSpikeFired=1*mV;beginRefrac=1*mV;inhibitionVoltage=prelimV',refractory=8*ms,clock=Clock(dt=self.generalClockDt))
				@network_operation(dt=self.generalClockDt)
				def additionToNetwork(): 
					neuronIndex = self.neuronIndex
					spikeIntervalCounter = (floor(self.time/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10
					def timePeriodAndRefractoryCalcs():
						self.time = self.time + self.timeStepInterval
						self.refractoryPointCounter = self.refractoryPointCounter + self.timeStepInterval	
						tT = (Decimal(format(self.time, '.1f')))
						sIU = (Decimal(format(spikeIntervalUnformatted, '.1f')))
						self.spikeIntervalCounter = int(tT/sIU)

						# if statement below corrects for offset of spiketimes starting at .1 sec.
						if self.time >= .1:
							self.neuronIndexCounter = self.neuronIndexCounter + self.timeStepInterval

						self.refractoryPointSwitch = False
						if self.refractoryPointCounter >= spikeIntervalUnformatted:
							self.refractoryPointCounter = 0.000
							self.refractoryPointSwitch = True
							self.lastSpikeInterval = self.time-self.timeStepInterval

						#print '__+++===','time',self.testTime,'refractoryPointCounter',self.testRefractoryPointCounter

					def WeightChangeCalculation(neuronIndex, spiketimes):
						for inputNeuronIndex in range(numberOfPixels):
							## General STDP learning rule implementation ##
							# Find SpikePreSyn if it exists
							# Note: could replace this with use of dictionary data structure for lookups if convenient
							# for processing time later

							# SpikePreSyn not found than make it max distance from SpikePostSyn
							# TODO SpikePreSyn = 0.1 # (Max spike time interval distance)
							# round up to nearest 100ms interval
							SpikePreSyn = self.time*second#(math.ceil(self.time*10)*.1)*second
							# SpikePostSyn not found than make it max distance from SpikePreSyn
							# TODO SpikePostSyn = 0.0
							#SpikePostSyn = 0*ms

							# .003 used below as estimation of value to use when no post syn spike is found and weights should
							# decrease.  Based on the formulas the weights would not decrease, as they are intended to, unless
							# a close post syn spike is used in a prior to the pre syn spike position
							SpikePostSyn = SpikePreSyn-(.001*second)#(math.floor(self.time*10)*.1)*second

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

								#print 'compare~~\t',CurrentSpikeNueron,inputNeuronIndex,CurrentSpikeTime,(self.time*second-.1*second),preSynSpikeFound

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
								CurrentSpikeTime = spikeCollection[1][i]#*10
								
								#print 'CurrentSpikeNueron',CurrentSpikeNueron,'neuronIndex',neuronIndex,'CurrentSpikeTime',CurrentSpikeTime,'self.time*second',self.time*second

								# exit loop once values below current time elapsed have all been checked
								# Disabled due to spikeCollection not being sorted and causing break too early
								#if CurrentSpikeTime > (self.time*second):
								#	break

								# (self.time*second-.1*second) is to check if in relevant time window below.
								# Note: that may not be a good cut off and I should check it
								# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
								# added (.1+.008)
								if CurrentSpikeNueron == neuronIndex and CurrentSpikeTime >= (self.time*second-(.1+.008)*second) and CurrentSpikeTime <= (self.time*second):
									SpikePostSyn = CurrentSpikeTime	
									postSynSpikeFound = True	

							if preSynSpikeFound == False or postSynSpikeFound == False:
								SpikePostSyn = SpikePreSyn-(.001*second)
								'''# Todo: watch more for spikes in a relevant time frame but are actually from different synapses and would
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
									#print 'F F'''

							#SpikePreSyn = 500*ms
							#SpikePostSyn = 700*ms
							
							#print 'neuronIndex', neuronIndex, 'inputNeuronIndex',inputNeuronIndex,'self.time', self.time, 'ADDS.v2',ADDS.v2, 'preSynSpikeFound',preSynSpikeFound,'postSynSpikeFound',postSynSpikeFound,'SpikePostSyn', SpikePostSyn, 'SpikePreSyn', SpikePreSyn

							# Find DeltaW
							DeltaW, DeltaSpikeTime = returnDeltaW(SpikePreSyn, SpikePostSyn)  

							# Find new weight
							#WOld = W[neuronIndex][inputNeuronIndex];
							WOld = dend[neuronIndex].w[inputNeuronIndex]

							# if statement below skips W change until initial spikes can be formed
							if self.time>.1: 
								NewW = returnNewW(WOld, DeltaW, DeltaSpikeTime)
								if NewW < 0: NewW = 0
								elif NewW > 1: NewW = 1
								#W[neuronIndex][inputNeuronIndex] = NewW
								dend[neuronIndex].w[inputNeuronIndex] = NewW*volt
								#print 'time',self.time,'NewW > 0',NewW,(NewW > 0),dend[neuronIndex].w[inputNeuronIndex]
								#print 'changes::','neuronIndex',neuronIndex,'inputNeuronIndex',inputNeuronIndex,'WOld',WOld,'dend[neuronIndex].w[inputNeuronIndex]',dend[neuronIndex].w[inputNeuronIndex],'DeltaW',DeltaW,'DeltaSpikeTime',DeltaSpikeTime

								#dend[neuronIndex].w[inputNeuronIndex] = W[neuronIndex][inputNeuronIndex]*volt
							else: 
								dend[neuronIndex].w[inputNeuronIndex] = WOld*volt
							#print 'DeltaW',DeltaW,'WOld',WOld,'WNew',W[neuronIndex][inputNeuronIndex]
							# Reuse existing values below instead of recalculating
							#somaDirect[neuronIndex].w[inputNeuronIndex] = dend[neuronIndex].w[inputNeuronIndex]

					def returnDeltaW(SpikePreSyn, SpikePostSyn):
						DeltaSpikeTime = SpikePreSyn - SpikePostSyn 
						DeltaSpikeTime = DeltaSpikeTime/second # remove second units
						# In thesis article it explains if DeltaT = 0 then DeltaW = 0
						DeltaW = 0
						# changed below line from that content in the article for a workaround/temp fix
						#if DeltaSpikeTime < 0:
						if DeltaSpikeTime < 0:
							DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus)))
							# for testing
							#DeltaW = 1
						elif DeltaSpikeTime > 0:
							DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus)))
							# for testing
							#DeltaW = -1

						return DeltaW, DeltaSpikeTime

					def returnNewW(WOld, DeltaW, DeltaSpikeTime):
						WOld = WOld/volt # remove volt unit
						WMin = 0;WMax = 1;
						# Check for inhibition occurence.  DeltaSpikeTime > 0 is used to represent inhibition synapse presence
						if DeltaSpikeTime > 0:
							WMin = -1;WMax = 0;	
						#print 'DeltaW', DeltaW
						## Relaxation rule implementation ##
						if DeltaW < 0:	
							WNew = WOld + ((LearningRate * (DeltaW * (WOld - WMin)))*self.negativeWeightReinforcement)
							
						elif DeltaW >= 0:				
							WNew = WOld + ((LearningRate * (DeltaW * (WMax - WOld)))*self.positiveWeightReinforcement)

						return WNew;
						##	

					def tauDCalc(neuronIndex, dendObject, W):
						# Weight loop
						for WIndex in range(len(W[0][:])):
							if abs(W[neuronIndex][WIndex]) <= 1:
								dendObject[neuronIndex].tau[WIndex] = (tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)) * ms# .001 is scaling factor found in dirac method

						#return dendObject[neuronIndex].tau

					def resistanceCalc(neuronIndex, dendObject, R):
						tauM = self.tauM
						#Resistance loop
						for RIndex in range(len(R[0][:])):				
							if (dendObject[neuronIndex].tau[RIndex]*.001 == tauM*second):
								# avoid a division by 0 issue
								tauM = tauM - .000001

							dendObject[neuronIndex].r[RIndex] = (((((dendObject[neuronIndex].tau[RIndex]*.001)/ms)*neuronFiringThreshold) / Rm) * ((tauM / ((dendObject[neuronIndex].tau[RIndex]/ms)*.001) ) ** (tauM / (tauM - ((dendObject[neuronIndex].tau[RIndex]/ms)*.001) ))))*volt

						#return dendObject[neuronIndex].r	

					def dendCalcs(neuronIndex, ADDSObj, dendObj, spiketimes, evaluationActive):
						# Below sequentially Dirac, Tau, then Resistance are calculated every end of a spike-time interval.
						# The resulting Dend I is added to the Um calc for the ADDS soma.

						# Dirac
						dendObj[neuronIndex].dirac = diracCalc(dendObj, neuronIndex, spiketimes)

						# Initialize weights
						if self.time == 0.000:
							dend[neuronIndex].w = W[neuronIndex]*volt
							#for inputNeuronIndex in range(numberOfPixels):
							#	dend[neuronIndex].w[inputNeuronIndex] = W[neuronIndex][inputNeuronIndex]*volt
							#print 'W::',W

						if (evaluationActive==False and self.refractoryPointSwitch==True):
							# Only change weights of neuron fired
							if ADDSObj.UmSpikeFired[neuronIndex] == 1*mV:
								# Weights
								WeightChangeCalculation(neuronIndex, spiketimes)	
								#print 'w',' time',self.time,'\t',dendObj[0].w,'\t',dendObj[1].w,'\t',dendObj[2].w,'\t',dendObj[3].w
							# Tau
							tauDCalc(neuronIndex, dendObj, self.W)
							# Resistance
							resistanceCalc(neuronIndex, dendObj, self.R)

						# TODO: check do I need additional loop below?
						for indexOfDend in range(dictionaryLongitude):
							ADDSObj.DendI[indexOfDend] = sum(dendObj[indexOfDend].v[:])*self.dendCalcScaling

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

					def diracCalc(IDend, neuronIndex, spiketimes):
						time = self.time
						e = math.e
						lastSpikeInterval = self.lastSpikeInterval
						dend = IDend
						# This method only calculates dirac
						# dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
						# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
						# though

						#print ':::range(len(IDend[neuronIndex][:])\t','self.time\t',self.time,'\tneuronIndex\t',neuronIndex,'\t',IDend,'\t',IDend[neuronIndex]

						dendGroup = [None]*len(dend[neuronIndex][:])
						for IdIndex in range(len(dend[neuronIndex][:])):
							# set tpresyn initially too far to trigger dirac
							tPreSyn = time + 1
							for presynInd in range(shape(spiketimes)[0]):
								comparedSpikeTime = spiketimes[presynInd][1]

								if comparedSpikeTime > (lastSpikeInterval + spikeIntervalUnformatted):
									break

								#print ':::',comparedSpikeTime,self.lastSpikeInterval,spikeIntervalUnformatted
								##print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'spiketimes[presynInd][0] == IdIndex',(spiketimes[presynInd][0] == IdIndex),'IdIndex',IdIndex,'comparedSpikeTime > self.lastSpikeInterval',(comparedSpikeTime > self.lastSpikeInterval),'comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted)',(comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted))
								'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
								# 0.0 to 300 ms, as the input is intended'''
								# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
								#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > lastSpikeInterval and comparedSpikeTime <= (lastSpikeInterval + spikeIntervalUnformatted):
								#print 'IdIndex',IdIndex,'time',time,'comparedSpikeTime',comparedSpikeTime,'time == comparedSpikeTime',(time == comparedSpikeTime),'(Decimal(format(time, .8f)) == Decimal(format(comparedSpikeTime, .8f)))',(Decimal(format(time, '.8f')) == Decimal(format(comparedSpikeTime, '.8f'))),Decimal(format(time, '.8f')),Decimal(format(comparedSpikeTime, '.8f'))
								#if comparedSpikeTime > .89:
								#	print 'compare',(math.floor(1000*Decimal(format(time, '.3f')))*.001 == math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001),math.floor(1000*Decimal(format(time, '.3f')))*.001,math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001
								if spiketimes[presynInd][0] == IdIndex and (math.floor(1000*Decimal(format(time, '.3f')))*.001 == math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001):
									## adding below if statement to avoid constant input stimulus in evaluation, that seem more accurate to how the article did it's evaulation.
									# normalize t to count just time within spike interval.  Formula from article just uses t in a way
									# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
									#tNorm = time - (floor((time/.001)*.01) * .1)
									#print 'tNorm',tNorm
									#if tNorm <= .01:# or evaluationActive == False:
									#if True:
										#tPreSyn = comparedSpikeTime
									tPreSyn = comparedSpikeTime
									#print '++++++++++',time,'+++++++++++','spiketimesnumber',spiketimes[presynInd][0],'neuronIndex',neuronIndex,'spiketimestime',(math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001),'comparedSpikeTime',comparedSpikeTime

							Dt = Decimal(format(time, '.8f')) - Decimal(format(tPreSyn, '.8f'))
							#print '++++++++Dt++++++++',Dt,'Decimal(format(time, .8f))',Decimal(format(time, '.8f')),'Decimal(format(tPreSyn, .8f))',Decimal(format(tPreSyn, '.8f'))
							# self.time == lines below are a workaround for initialization values
							if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
								Dt = 1.0
							##print 'neuronIndex',neuronIndex,'tPreSyn000:',tPreSyn,'Dt',Dt
							
							# simplify dirac until later.  TODO: try more comple dirac
							#if (t > -(Dt/2) and t < (Dt/2)):
							# Seems to me what is looked for here is that the post synapse (output neurons) is after the pre synapse (input neurons)
							if Dt >= 0.0:
								if self.time == 0.121 or self.time == 0.421 or self.time == 0.721 or self.time == 1.021:
									DiracFun = 1.00
								if Dt != 0: 
									# A closer post to pre synapse creates a larger dirac signal for greater weight I'd suppose than
									# a further distance below in DiracFun
									#Dt = Dt*1000
									DiracFun = 1/Dt
								else:
									DiracFun = 1#1000
								DiracFun = 1
								#dendGroup[IdIndex] = float(DiracFun)*volt*.001*self.diracScaling #add .001 as scaling factor.  
								# Scaling factor of .001 could account for 1ms being present instead of .001 as it is computed, difference of
								# 1/1 compared to 1/.001
								dendGroup[IdIndex] = float(DiracFun)*volt*self.diracScaling
							else:
								DiracFun = 0
								dendGroup[IdIndex] = float(DiracFun)*volt*.001*self.diracScaling						

						return dendGroup										

					def checkForResets(neuronIndex, dendObj, somaDirectObj):
						if self.v2[neuronIndex] == 10*mV:
							self.v2[neuronIndex] = 0*mV	
							'''dendObj[neuronIndex].dirac = [0*volt]*len(dendObj[neuronIndex][:])
							dendObj[neuronIndex].v = [0*volt]*len(dendObj[neuronIndex][:])	'''						
							# Activate lateral inhibition upon a spike to inhibit other neurons from spiking.  As brian's 
							# Inhibition inhibits before dend input can increase the signal to balance things out.
							#ADDS.lateralInhibition()
						elif self.v[neuronIndex] != -0.002 * mV:
							self.v2[neuronIndex] = self.v[neuronIndex]

						if somaDirectObj.v2[neuronIndex] == 10*mV:
							somaDirectObj.v2[neuronIndex] = 0*mV	

						if (self.beginRefrac[neuronIndex] == 1*mV) and (self.relRefracTimeActivated[neuronIndex] == False):
							self.beginRefrac[neuronIndex] = 0*mV
							self.absRefracTime[neuronIndex] = self.time
							self.relRefracTime[neuronIndex] = self.time
							self.absRefracTimeActivated[neuronIndex] = True
							self.relRefracTimeActivated[neuronIndex] = True
						# Absolute refrac
						if self.absRefracTimeActivated[neuronIndex] == True:
							'''testADDS.v[neuronIndex] = -0.002*mV
							testADDS.DendI[neuronIndex] = -0.002*mV
							testADDS.SynI[neuronIndex] = -0.002*mV
							testDend[neuronIndex].dirac = [0*volt]*len(dend[neuronIndex][:])'''
							if self.time >= (self.absRefracTime[neuronIndex] + self.absRefracTimeDuration):
								self.absRefracTimeActivated[neuronIndex] = False
						#print 'neuronIndex', neuronIndex, 'abs',self.testAbsRefracTimeActivated[neuronIndex],'self.testTime',self.testTime,'self.beginRefrac[neuronIndex]',self.beginRefrac[neuronIndex],'self.relRefracTimeDuration/ms',(self.absRefracTimeDuration/ms),'self.testAbsRefracTime[neuronIndex] + self.absRefracTimeDuration',(self.testAbsRefracTime[neuronIndex] + self.absRefracTimeDuration)
						# Relative refrac
						if self.relRefracTimeActivated[neuronIndex] == True:
							'''if testADDS.v[neuronIndex] > -0.002 * mV:
								testADDS.v[neuronIndex] = -0.002*mV
								testADDS.DendI[neuronIndex] = -0.002*mV
								testADDS.SynI[neuronIndex] = -0.002*mV
								testDend[neuronIndex].dirac = [0*volt]*len(dend[neuronIndex][:])'''
							if self.time >= (self.relRefracTime[neuronIndex] + self.relRefracTimeDuration):
								self.relRefracTimeActivated[neuronIndex] = False
						#print 'neuronIndex', neuronIndex, 'rel',self.testRelRefracTimeActivated[neuronIndex],'self.testTime',self.testTime,'self.testRelRefracTime[neuronIndex]',self.testRelRefracTime[neuronIndex],'self.relRefracTimeDuration/ms',(self.relRefracTimeDuration/ms),'self.testRelRefracTime[neuronIndex] + self.relRefracTimeDuration',(self.testRelRefracTime[neuronIndex] + self.relRefracTimeDuration)

					def printOutputForTesting():
						print 'neuronIndex',self.neuronIndex,'time',self.time
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
						print '3',dend[3].r[:]			

					def sortVersionOutForTesting(ADDSObj, dendObj):
						print 'ADDSObj.v','time',self.time,ADDSObj.v[0],ADDSObj.v[1],ADDSObj.v[2],ADDSObj.v[3]
						print 'sum(ADDSObj.DendI)',(ADDSObj.DendI[0]),(ADDSObj.DendI[1]),(ADDSObj.DendI[2]),(ADDSObj.DendI[3])
						print 'sum(ADDSObj.SynI)',sum(ADDSObj.SynI[0]),sum(ADDSObj.SynI[1]),sum(ADDSObj.SynI[2]),sum(ADDSObj.SynI[3])						

					def printEvalOutputForTesting(ADDSObj, dendObj):
						'''print 'dendObj.dirac','time',self.time,'\t',dendObj[0].dirac,dendObj[1].dirac,dendObj[2].dirac,dendObj[3].dirac
						print 'dendObj.v',dendObj[0].v,dendObj[1].v,dendObj[2].v,dendObj[3].v
						print 'sum(dendObj[neuronIndex])',sum(dendObj[0].v[:]),sum(dendObj[1].v[:]),sum(dendObj[2].v[:]),sum(dendObj[3].v[:])							
						print 'testSomaDirect.v',testSomaDirect.v[0],testSomaDirect.v[1],testSomaDirect.v[2],testSomaDirect.v[3]
						print 'sum(testSomaDirect[neuronIndex])',sum(testSomaDirect.v[0]),sum(testSomaDirect.v[1]),sum(testSomaDirect.v[2]),sum(testSomaDirect.v[3])							
						print 'testSomaDirect.summedWandDirac',testSomaDirect.summedWandDirac[0],testSomaDirect.summedWandDirac[1],testSomaDirect.summedWandDirac[2],testSomaDirect.summedWandDirac[3]					'''
						print 'ADDSObj.v',ADDSObj.v[0],ADDSObj.v[1],ADDSObj.v[2],ADDSObj.v[3]
						print 'dendObj.dirac','time',self.time,'\t',dendObj[0].dirac,dendObj[1].dirac,dendObj[2].dirac,dendObj[3].dirac
						print 'dendObj.v',dendObj[0].v,dendObj[1].v,dendObj[2].v,dendObj[3].v
						print 'v',dendObj[0].v,'\t',dendObj[1].v
						print 'r',dendObj[0].r,'\t',dendObj[1].r
						print 'w',' time',self.time,'\t',dendObj[0].w,'\t',dendObj[1].w,'\t',dendObj[2].w,'\t',dendObj[3].w
						print 'dirac',dendObj[0].dirac,'t',dendObj[1].dirac
						print 'tau',dendObj[0].tau,'t',dendObj[1].tau
						print 'units',(-dendObj[1].v[0]/mV),((dendObj[1].r[0])/volt),((dendObj[1].w[0])/volt),(dendObj[1].dirac[0]/volt),(dendObj[1].tau[0]/second)
						dendObj[1].vTest[0] = (((-dendObj[1].v[0]/mV)+(((dendObj[1].r[0])/volt)*((dendObj[1].w[0])/volt)*(dendObj[1].dirac[0]/volt)))/(dendObj[1].tau[0]/second))*volt#((-dendObj[1].v[0]/volt)+((dendObj[1].r[0])/volt))/(dendObj[1].tau[0]/ms)*volt
						dendObj[1].vTest2[0] = -dendObj[1].v[0]
						dendObj[1].vTest3[0] = ((dendObj[1].r[0]/volt)*(dendObj[1].w[0]/volt)*(dendObj[1].dirac[0]))/(dendObj[1].tau[0]/second)
						print 'vTest',dendObj[1].vTest[0]
						print 'vTest2',dendObj[2].vTest2[0],dendObj[1].v[0],dendObj[1].v[1],dendObj[1].v[2],dendObj[1].v[3],dendObj[1].v[4]
						print 'vTest3',dendObj[3].vTest3[0]
						#print 'rs equ',(-dendObj[0].v+((dendObj[0].r)*(dendObj[0].w)*dendObj[0].dirac)/(dendObj[0].tau)),'\t',(-dendObj[1].v+((dendObj[1].r)*(dendObj[1].w)*dendObj[1].dirac)/(dendObj[1].tau))
						print 'sum(ADDSObj.DendI)',(ADDSObj.DendI[0]),(ADDSObj.DendI[1]),(ADDSObj.DendI[2]),(ADDSObj.DendI[3])
						print 'sum(ADDSObj.SynI)',sum(ADDSObj.SynI[0]),sum(ADDSObj.SynI[1]),sum(ADDSObj.SynI[2]),sum(ADDSObj.SynI[3])
						print 'ADDSObj.UmSpikeFired',ADDSObj.UmSpikeFired[0],ADDSObj.UmSpikeFired[1],ADDSObj.UmSpikeFired[2],ADDSObj.UmSpikeFired[3]
						#print 'Prelim Weights\n0',dendObj[0].w[:]
						#print '1',dendObj[1].w[:]
						#print '2',dendObj[2].w[:]
						#print '3',dendObj[3].w[:]
						#print 'Prelim Tau\n0',dendObj[0].tau[:]
						#print '1',dendObj[1].tau[:]
						#print '2',dendObj[2].tau[:]
						#print '3',dendObj[3].tau[:]
						#print 'Prelim Res\n0',dendObj[0].r[:]
						#print '1',dendObj[1].r[:]
						#print '2',dendObj[2].r[:]
						#print '3',dendObj[3].r[:]
						print '____________________'						

					def evaluateClassifier(ADDSObj):
						# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
						# (self.testSpikeIntervalCounter-1) is to correct for refractoryPointSwitch occuring after spikeInterval it addresses.
						#print 'refractoryPointSwitch = ', self.refractoryPointSwitch
						#print 'refractoryPointSwitch =', self.testRefractoryPointSwitch
						#if self.refractoryPointSwitch == true and (self.testSpikeIntervalCounter > 0):
						if self.refractoryPointSwitch == True:# and (self.testSpikeIntervalCounter > 0):
							'''print 'self.testTime',self.testTime,'refrac on'
							tT = (Decimal(format(self.testTime, '.1f')))
							sIU = (Decimal(format(spikeIntervalUnformatted, '.1f')))
							timeDiff = tT/sIU
							print 'recalc testT',(floor(self.testTime/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10,'floor(self.testTime/spikeIntervalUnformatted)',floor(self.testTime/spikeIntervalUnformatted),'floor with decimal',int(timeDiff)'''
						#if self.testSpikesFiredInInterval[neuronIndex][self.testSpikeIntervalCounter-1] == False:
							# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
							# the +1 in (self.timeotalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
							# * A change was made to the scoring to not count the 0.0-0.1 second period because the input spike generator does not start until .1
							# seconds and the first occurences of output spikes should be monitored looking at .2 seconds to see if any occured in seconds .1-.2 .
							if (self.spikeIntervalCounter >= 2) and (self.spikeIntervalCounter <= (self.totalSpikeIntervals+1)):
								for neuronIndex in range(dictionaryLongitude):				
									#print 'eval','self.testSpikeIntervalCounter\t',self.testSpikeIntervalCounter
									if (ADDSObj.UmSpikeFired[neuronIndex] == 1*mV):# and (self.testSpikeIntervalCounter < self.totalSpikeIntervals):
										if (self.correctSpikes[neuronIndex][(self.spikeIntervalCounter-1)] == 1):
											self.truePositiveSpikeResults = self.truePositiveSpikeResults + 1	
											print 'TP found\t','self.testSpikeIntervalCounter-1\t',self.spikeIntervalCounter-1,'neuronIndex\t',neuronIndex
										else:
											self.falsePositiveSpikeResults = self.falsePositiveSpikeResults + 1	
										#self.testSpikesFiredInInterval[neuronIndex][self.testSpikeIntervalCounter] = True	
									elif (ADDSObj.UmSpikeFired[neuronIndex] == 0*mV):# and (self.testSpikeIntervalCounter < self.totalSpikeIntervals):
										if (self.correctSpikes[neuronIndex][(self.spikeIntervalCounter-1)] == 1):
											self.falseNegativeSpikeResults = self.falseNegativeSpikeResults + 1		
										else:
											self.trueNegativeSpikeResults = self.trueNegativeSpikeResults + 1	
							for neuronIndex in range(dictionaryLongitude):											
								ADDSObj.UmSpikeFired[neuronIndex] = 0*mV
									#print 'results',self.truePositiveSpikeResults,self.falsePositiveSpikeResults,self.trueNegativeSpikeResults,self.falseNegativeSpikeResults

					def intitializeTrainedModelParameters(dendObj):
						# Below values that have been created through prior training are used for w, tau, and r
						testW = [[0.90433194, 0.6139531, 0.50387484, 0.55220372, 0.51213536, 0.85443374, 0.99955922, 0.5039825, 0.73091913, 0.9780236, 0.5241028, 0.71571812, 0.93782861, 0.51210244, 0.73074697],
						[0., 0., 0.03412608, 0., 0.90455366, 0.78683668, 0., 0.95912629, 0.7282637, 0., 0.78548583, 0.78935491, 0.03193823, 0.00609877, 0.17287094],
						[0.4474444, 0., 0.98135641, 0., 0.96315942, 0., 0., 0., 0.15930208, 0., 0.77299245, 0., 0., 0.71739497, 0.02804206],
						[0., 0., 0.99815102, 0., 0.9239562, 0., 0., 0.32862838, 0.29682383, 0., 0.85108903, 0., 0., 0., 0.6687179]]

						testTauD = [[4.67870579, 12.80931328, 15.89150457, 14.53829579, 15.66020997, 6.07585538, 2.0123418, 15.88849011, 9.53426444, 2.6153392, 15.32512158, 9.95989264, 3.74079895, 15.66113156, 9.53908472],
						[9.38715598, 12.1194624, 9.34555739, 6.3921125, 5.61075905, 4.06988445, 9.88122329, 3.54626289, 6.41301412, 7.44992802, 10.1151254, 3.96914614, 3.41772316, 4.7844492, 2.40562813],
						[14.7112847, 14.11913186, 2.72683686, 10.43654572, 3.43626265, 14.20735373, 9.15843206, 4.26396028, 3.39998086, 7.5902114, 10.85008977, 7.25459199, 9.86718901, 3.20308826, 10.16608643],
						[15.97156807, 4.01171718, 2.07503494, 9.24935852, 5.08599601, 13.72643044, 13.71622469, 7.36132416, 5.76486557, 7.2621342, 8.04307852, 11.68604314, 5.89535223, 10.07822174, 15.44403097]]

						testR = [[5.28618924, 7.07014396, 7.67118426, 7.41054119, 7.62696131, 5.62546396, 4.55408732, 7.67060873, 6.39667963, 4.73398033, 7.56265998, 6.48674751, 5.04468982, 7.62713777, 6.39770456],
						[6.36534465, 6.93170042, 6.35646422, 5.69956225, 5.51481581, 5.13093667, 6.47016506, 4.99283791, 5.70442879, 5.94150634, 6.51938451, 5.10472189, 4.95819508, 5.31264187, 4.67271447],
						[7.44412582, 7.32882233, 4.76605631, 6.58661617, 4.96321112, 7.34606257, 6.31640753, 5.1810014, 4.95338848, 5.97297472, 6.67245869, 5.8974608, 6.46720377, 4.89962516, 6.53007485],
						[7.68646215, 5.11581952, 4.57335711, 6.33589406, 5.38732016, 7.25180984, 7.2498024, 5.92156072, 5.55171007, 5.89916649, 6.07367171, 6.84386832, 5.58276752, 6.51163569, 7.58550996]]
						#print 't:',self.time,'neuronIndex',neuronIndex,'main calc resistanceCalc()', testR

						for indexOfDend in range(dictionaryLongitude):
							# TODO see if some of these do not need to be recomputed every time
							dendObj[indexOfDend].w = testW[indexOfDend]*volt*self.weightScaling
							dendObj[indexOfDend].tau = testTauD[indexOfDend]*ms # unclear if dividing by tau in ms is e.x. /.02 or /20 but it is assumed to be 20, therefore no ms conversion here
							dendObj[indexOfDend].r = testR[indexOfDend]*mV	# volt unit is cancelled out in the equation anyhow, doesn't matter if it is volt or mV due to being cancelled.  Having mv could cause *.001 that is now wanted

					def lateralInhibition(ADDSObj):
						'''spikeOccured = False;inhibitionVoltage = None
						for neuronIndex3 in range(dictionaryLongitude):
							if ADDSObj.UmSpikeFired[neuronIndex3] == 1*mV: 
								spikeOccured = True
								inhibitionVoltage = ADDSObj.inhibitionVoltage[neuronIndex3]#*self.lateralInhibitScaling'''

						tNorm = self.time - (floor((self.time/.001)*.01) * .1)
						preliminaryV = [None]*dictionaryLongitude
						greatestMembraneVoltage = 0.0
						greatestNeuronNumber = None
						neuronNumberFired = None
						neuronNumberAboutToFire = None
						
						for neuronIndex in range(dictionaryLongitude):
							preliminaryV[neuronIndex] = ADDSObj.prelimV[neuronIndex]
						for neuronIndex in range(dictionaryLongitude):
							for neuronIndex2 in range(len(preliminaryV)):
								if (neuronIndex != neuronIndex2) and (preliminaryV[neuronIndex2] >= greatestMembraneVoltage): 
									#print 'comparison','neuronIndex',neuronIndex,' ',preliminaryV[neuronIndex],'neuronIndex2',neuronIndex2,' ',preliminaryV[neuronIndex2],'greatestMembraneVoltage',greatestMembraneVoltage,preliminaryV[neuronIndex2] >= greatestMembraneVoltage
									greatestMembraneVoltage = preliminaryV[neuronIndex2]
									greatestNeuronNumber = neuronIndex
								if ADDSObj.UmSpikeFired[neuronIndex2] == 1*mV: neuronNumberFired = neuronIndex2

							if (greatestMembraneVoltage >= (10*mV)) and (neuronNumberFired == None):
								self.voltageForInhibition = greatestMembraneVoltage - (10*mV)
								neuronNumberAboutToFire = greatestNeuronNumber
								self.inhibitionReduction = 1.0
							elif (greatestMembraneVoltage < (10*mV)) and (neuronNumberFired == None): 
								self.voltageForInhibition = greatestMembraneVoltage
							else:
								self.voltageForInhibition = 70*mV#(70*mV/self.inhibitionReduction)
								if self.inhibitionReduction < 150: self.inhibitionReduction += 0.01
							#print 'time',self.time,'neuronIndex',neuronIndex,'neuronIndex2',neuronIndex2,'self.voltageForInhibition',self.voltageForInhibition,(70*mV/self.inhibitionReduction),'self.inhibitionReduction',self.inhibitionReduction

						for neuronIndex in range(dictionaryLongitude):		
							if neuronIndex != neuronNumberAboutToFire and neuronIndex != neuronNumberFired:
								#if self.voltageForInhibition > 0: preliminaryV[neuronIndex] -= self.voltageForInhibition
								#else: preliminaryV[neuronIndex] += self.voltageForInhibition
								preliminaryV[neuronIndex] -= abs(self.voltageForInhibition)
								ADDSObj.v[neuronIndex] = preliminaryV[neuronIndex]

							'''for neuronIndex2 in range(dictionaryLongitude):
								#inhibitionVoltage = ADDSObj.prelimV[neuronIndex2]*self.lateralInhibitScaling
								#if ADDSObj.UmSpikeFired[neuronIndex2] == 1*mV: 
								#	inhibitionVoltage = ADDSObj.inhibitionVoltage[neuronIndex2]#*self.lateralInhibitScaling								
								if neuronIndex != neuronIndex2 and self.voltageForInhibition > 0:
									preliminaryV[neuronIndex] -= inhibitionVoltage'''
							#print 'self.time',self.time,'neuronIndex',neuronIndex,'self.voltageForInhibition',self.voltageForInhibition,'ADDSObj.prelimV[neuronIndex]',ADDSObj.prelimV[neuronIndex],'preliminaryV[neuronIndex]',preliminaryV[neuronIndex],'ADDSObj.v[neuronIndex]',ADDSObj.v[neuronIndex],'ADDSObj.UmSpikeFired[neuronIndex]',ADDSObj.UmSpikeFired[neuronIndex]
							#if ADDSObj.UmSpikeFired[neuronIndex] != 1*mV: ADDSObj.v[neuronIndex] = preliminaryV[neuronIndex]

					def mainSimulationCalcs(ADDSObj, dendObj, somaDirectObj, spiketimes, evaluationActive):
						# dend then somaDirect calcs are done which are then used to set lat inhib.
						# Soma Um calcs are done automatically using equations entered for brian
						# once dend and somaDirect are updated
						tNorm = self.time - (floor((self.time/.001)*.01) * .1)
						#print 'tNorm',tNorm
						
						timePeriodAndRefractoryCalcs()

						#print 'time',self.time,'brian time',self.t

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

							lateralInhibition(ADDSObj)

							for neuronIndex in range(dictionaryLongitude): checkForResets(neuronIndex, dendObj, somaDirectObj)

							evaluateClassifier(ADDSObj)

							'''if tNorm <= 0.003 or tNorm >= .099:
								printEvalOutputForTesting(ADDS, dend)
							else:
								sortVersionOutForTesting(ADDS, dend)'''

					if self.evaluateClassifier == False:
						mainSimulationCalcs(ADDS, dend, somaDirect, self.trainingSpikeTimes, False)
					else:
						mainSimulationCalcs(ADDS, dend, somaDirect, self.spiketimes, True)					

				self.contained_objects.append(additionToNetwork)	

			def OutputEvaluationResults(self, dendObj):
				print 'Final Weights\n',dendObj[0].w[:]/volt
				print dendObj[1].w[:]/volt
				print dendObj[2].w[:]/volt
				print dendObj[3].w[:]/volt
				print 'Final Tau\n',dendObj[0].tau[:]/ms
				print dendObj[1].tau[:]/ms
				print dendObj[2].tau[:]/ms
				print dendObj[3].tau[:]/ms
				print 'Final Res\n',dendObj[0].r[:]/mV
				print dendObj[1].r[:]/mV
				print dendObj[2].r[:]/mV
				print dendObj[3].r[:]/mV		
				print '\n'
				print '+++ Results +++'
				print 'Spike results: TP:\t',self.truePositiveSpikeResults,'\tFP:\t',self.falsePositiveSpikeResults,'\tTN:\t',self.trueNegativeSpikeResults,'\tFN:\t',self.falseNegativeSpikeResults
				print 'totalSpikeIntervalsTested:\t',self.totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
				print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(self.truePositiveSpikeResults, '.1f'))/Decimal(format(self.totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
				print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',(Decimal(format(self.truePositiveSpikeResults, '.1f'))+Decimal(format(self.trueNegativeSpikeResults, '.1f')))/(Decimal(format(self.totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
				print '+++++++++++++++'				

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

		ADDS.OutputEvaluationResults(dend)

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
