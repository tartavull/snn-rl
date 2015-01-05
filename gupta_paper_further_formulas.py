from architecture_further_formulas import *

class gupta_paper:
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	#print 'st info:\t',spiketimes[0:310,0],'\t__\t',spiketimes[0:310,1]
	testSpiketimes = spiketimes
	LIK = SpikeGeneratorGroup(N=15, indices=spiketimes[:,0], times=spiketimes[:,1]*ms)
	neuronIndex = 0
	t = 0.0000
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
	lastSpikeInterval = 0.0
	currentEpoch = 0
	refractoryPointCounter = 0.0
	#refractoryPointCounter = -0.008 # for workaround adjusting this variable
	#refractoryPointCounter = -0.03 # for workaround adjusting this variable
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

	print 'initial Weights\n',W

	def run_model(self):
		
		dictionary = self.dictionary
		spiketimes = self.spiketimes
		#tauMZ*dv/dt = -v + RmZ*(IdVZ) : volt
		#dv/dt = (-v + (RmZ/mV)*IdVZ)/(tauMZ) : volt
		eqs = Equations('''
			dv/dt = (-v+(IdVZ))/(20*ms) : volt (unless refractory)
			RmZ = 80*mV : volt
			tauMZ = 20*ms : second
			IdVZ : volt
	        V : volt
	        DendI : volt
	        SynI : volt
	        v2 : volt
	        dendriV1 : volt
	        dendriV2 : volt
	        dendriV3 : volt
	        dendriV4 : volt
	        dendriV5 : volt
	        dendriV6 : volt
	        dendriV7 : volt
	        dendriV8 : volt
	        dendriV9 : volt
	        dendriV10 : volt
	        dendriV11 : volt
	        dendriV12 : volt
	        dendriV13 : volt
	        dendriV14 : volt			
		    ''')			

		#			dv/dt = (-v + (IdVZ))/(tauMZ) : volt
		myclock=Clock(dt=1*ms)

		ADDS = NeuronGroup(N=4, model=eqs,threshold='v>.01*mV', reset='v=-0.002 * mV; dv=0; IdVZ = 0*mV;v2=10*mV',refractory=10*ms,clock=myclock)

		def returnUm(self):
			# Main logic of the neuron simulation model is performed here
			# Modules are performed below resulting in voltage for the neurons being calculated

			dictionary = self.dictionary
			neuronIndex = self.neuronIndex
			t = self.t
			tNorm = t - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
			spiketimes = self.spiketimes
			
			# Calculate tauD
			tauD = self.tauD
			tauD = tauDCalc(neuronIndex, tauD, W)
			#print 'main calc tauDCalc()', tauD

			# Calculate resistance
			R = self.R
			R = resistanceCalc(neuronIndex, tauD, R)
			#print 'main calc resistanceCalc()', R

			#Dedritic total post-synaptic current
			#Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(self.neuronIndex)
			Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, Id, self.IdSpikeFired, self.IdRefractoryPeriod, R, W, spiketimes)
			#print 'neuronIndex',neuronIndex,'self.t',self.t,'ADDS.t',ADDS.t,'Id[neuronIndex]',Id[neuronIndex]
			# Direct to soma 
			#Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, Id, self.IdSpikeFired, self.IdRefractoryPeriod, R, W, spiketimes)

			### Soma membrane potential ###
			Um[neuronIndex] = totalSomaMembranePotential(neuronIndex, Um, Id, Is, tNorm)

			# Refractory and spike evaluation: Is
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])			

			# Refractory and spike evaluation: Um
			Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])
			Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])

			timePeriodAndRefractoryCalcs()

			return Um[neuronIndex] * mV;	

		def refractoryPeriodEvaluation(Voltage, SpikeFired, RefractoryPeriod):
			# Check for spike
			# when identifying a spike, refractoryPeriodEndPoint is used as an alternative condition to refractoryPeriod being
			# over to allow a spike right at that end point

			# Refractory period
			#refractoryPeriodEndPoint = (self.t - 0.001) % spikeIntervalUnformatted <= 0.00001			
			if SpikeFired == True and self.refractoryPointSwitch == False:
				RefractoryPeriod = True
				#print 'refrac activated!'
			elif RefractoryPeriod == True and self.refractoryPointSwitch == True:
				# 0.001 added above for time step compatibility.   0.00001 instead of 0.0 used for some kind of 
				# modulo computational result offset in python which was producing a value from the modulo calc
				# just slightly over 0
				RefractoryPeriod = False
				#print 'refrac over!'
			#print 'self.refractoryPeriod == True and refractoryPeriodEndPoint', self.refractoryPeriod == True, refractoryPeriodEndPoint, (self.refractoryPeriod == True and refractoryPeriodEndPoint), self.t % (spikeIntervalUnformatted + 0.001) <= 0.00001
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

		def dentritePostSynapticCurrent(neuronIndex, IDend, IDendSpikes, IDendRefract, R, W, spiketimes):
			# Solving for Idend in the formula in the article yeilded the below equation
			t = self.t
			e = math.e

			for IdIndex in range(len(IDend[neuronIndex][:])):
				tauDen = tauD[neuronIndex][IdIndex]
				r = R[neuronIndex][IdIndex]
				w = W[neuronIndex][IdIndex]

				# set tpresyn initially too far to trigger dirac
				#tPreSyn = -t - 1
				tPreSyn = t + 1
				for presynInd in range(shape(spiketimes)[0]):
					comparedSpikeTime = spiketimes[presynInd][1]
					#print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval
					'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
					# 0.0 to 300 ms, as the input is intended
					# changed <= to < below to see if that works better
					if spiketimes[presynInd][0] == IdIndex and (comparedSpikeTime-spikeIntervalUnformatted) < t:'''
					# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
					#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > (self.lastSpikeInterval - spikeIntervalUnformatted):
					#print 'self.t',self.t,'neuronIndex',neuronIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0],'comparedSpikeTime',comparedSpikeTime#,'self.lastSpikeInterval',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'IdIndex',IdIndex,'spiketimes[presynInd][0]',spiketimes[presynInd][0]
					if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime >= (self.lastSpikeInterval - spikeIntervalUnformatted) and comparedSpikeTime <= self.lastSpikeInterval:
						tPreSyn = comparedSpikeTime

					# this is true? -> commented out logic below due to spiketimes perhaps not being in sorted order, a relevant entry could occur after an irrelevant one
					#if (comparedSpikeTime-spikeIntervalUnformatted) > t:
					# spike times are in order of earliest to latest times grouped with input pixel number.  therefore the below cutoff in the loop is fine.
					if comparedSpikeTime > self.lastSpikeInterval:
						break

				#print 'tPreSyn000:',tPreSyn
				Id2 = IDend[neuronIndex][IdIndex]
				Dt = t - tPreSyn
				# self.t == lines below are a workaround for initialization values
				if self.t == 0.121 or self.t == 0.421 or self.t == 0.721 or self.t == 1.021:
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
				# for testing .025 below is scaling factor to get spikeCoeff to be at least .011 which is enough to cause spikes.
				#SpikeModCoeff = (r*(w*.025))
				#SpikeModCoeff = (r*(w*.012))
				#SpikeModCoeff = (r*(w*.0024))
				#r = .00375
				#r = .0025
				#r = .0045
				SpikeModCoeff = (r*w)
				
				#print '***ni:\t',neuronIndex,'\tIdIndex\t',IdIndex,'\t***w:\t',w

				# correct for scaling
				tauDen = tauDen * .001
				
				# normalize t to count just time within spike interval.  Formula from article just uses t in a way
				# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
				tNorm = t - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
				#if -Dt/2<(tNorm or t?)<Dt/2:
				#if Dt <= 0.0:
				if Dt >= 0.0:
					#DiracFun = 1
					#DiracFun = 0.95
					#DiracFun = 1.00
					if self.t == 0.121 or self.t == 0.421 or self.t == 0.721 or self.t == 1.021:
						DiracFun = 1.00
					#SpikeModCoeff = (r*w*DiracFun)
					#SpikeModCoeff = .011
					#SpikeModCoeff = .005
					#SpikeModCoeff = .011
					#DiracFun = 0.66
					SpikeModCoeff = (SpikeModCoeff*DiracFun)
					#tauDen = .03	
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
				else:
					DiracFun = 0
					#SpikeModCoeff = (r*w*DiracFun)
					#SpikeModCoeff = .011
					#SpikeModCoeff = .005
					#SpikeModCoeff = .011
					SpikeModCoeff = (SpikeModCoeff*DiracFun)					
					#tauDen = .03										
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
				
				#print 'part1: \t',IDend
				#print 'factors\tSpikeModCoeff\t',SpikeModCoeff,'self.t',self.t,'neuronIndex',neuronIndex,'\tId2\t',Id2,'\ttPreSyn\t',tPreSyn,'\tDiracFun\t',DiracFun,'Dt <= 0.0',(Dt <= 0.0),'Dt',Dt,'\tneuronIndex\t',neuronIndex,'\tIdIndex\t',IdIndex,'-(SpikeModCoeff - Id2)',-(SpikeModCoeff - Id2),'-(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen))',-(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen))

				# Refractory and spike evaluation: Id
				#refractResults = refractoryPeriodEvaluation(IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex])
				#print 'part2: \t',refractResults
				#spikeEvalResults = spikeFiredEvaluation(refractResults[0], refractResults[1], refractResults[2])
				#print 'part3: \t',spikeEvalResults
				spikeEvalResults = IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex]
				if self.t > (self.lastSpikeInterval+0.02):
					refractResults = refractoryPeriodEvaluation(IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex])
					spikeEvalResults = spikeFiredEvaluation(refractResults[0], refractResults[1], refractResults[2])
				
				IDend[neuronIndex][IdIndex] = spikeEvalResults[0]
				IDendSpikes[neuronIndex][IdIndex] = spikeEvalResults[1]
				IDendRefract[neuronIndex][IdIndex] = spikeEvalResults[2]

				#print 'part1_2: \t', IDend
				
			return [IDend[neuronIndex], IDendSpikes[neuronIndex], IDendRefract[neuronIndex]]

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

		def totalSomaMembranePotential(neuronIndex, UMemb, IDend, ISoma, tNorm):
			### Soma membrane potential ###
			# Solving for Um in the formula in the article yeilded the below equation
			Um2 = UMemb[neuronIndex]

			SummedDendriteGroup = sum(IDend[neuronIndex])
			#print 'neuronIndex', neuronIndex, 'self.t', self.t, 'tNorm',tNorm,'SummedDendriteGroup', SummedDendriteGroup
			ADDS.IdVZ[neuronIndex] = SummedDendriteGroup*volt*10#*10 added 12/29/14 for scaling adjustment test # * 1 # the * 1000 is for a scaling adjustment
			SynapseToSoma = ISoma[neuronIndex]
			#print 'self.epochIndex\t',self.epochIndex,'\tneuronIndex\t',neuronIndex,'\tSummedDendriteGroup\t',SummedDendriteGroup,'\tIDend\t',IDend

			#self.SummedDendriteGroup = SummedDendriteGroup;
			#self.SynapseToSoma = SynapseToSoma;

			#NeuronInputCoeff = Rm*(SummedDendriteGroup)
			NeuronInputCoeff = SummedDendriteGroup

			newTSMP = -(NeuronInputCoeff - Um2) * (e ** (-tNorm/tauM)) + NeuronInputCoeff			

			#print 'self.t',self.t,'\tSummedDendriteGroup\t',SummedDendriteGroup,'neuronIndex\t',neuronIndex,'\tnewTSMP\t',newTSMP,'\ttestId[neuronIndex]\t',testId[neuronIndex]

			return newTSMP

		def WeightChangeCalculation():
			for inputNeuronIndex in range(numberOfPixels):
				## General STDP learning rule implementation ##
				# Find SpikePreSyn if it exists
				# Note: could replace this with use of dictionary data structure for lookups if convenient
				# for processing time later

				# SpikePreSyn not found than make it max distance from SpikePostSyn
				# TODO SpikePreSyn = 0.1 # (Max spike time interval distance)
				# round up to nearest 100ms interval
				SpikePreSyn = (math.ceil(self.t*10)*.1)*second
				# SpikePostSyn not found than make it max distance from SpikePreSyn
				# TODO SpikePostSyn = 0.0
				#SpikePostSyn = 0*ms
				SpikePostSyn = (math.floor(self.t*10)*.1)*second

				preSynSpikeFound = False
				postSynSpikeFound = False

				spikeCollection = spiketimes
				NumberOfSpikes = shape(spikeCollection)[0]
				for i in range(NumberOfSpikes):
					CurrentSpikeNueron = spikeCollection[i][0]
					CurrentSpikeTime = spikeCollection[i][1]*second

					# exit loop once values below current time elapsed have all been checked
					if CurrentSpikeTime > (self.t*second):
						print 'eeww1:\tspiketimes\t',spiketimes,'\tqq\t',CurrentSpikeTime
						break

					# (self.t*second-.1*second) is to check if in relevant time window below.
					# Note: that may not be a good cut off and I should check it
					# added < self.t just for testing
					if CurrentSpikeNueron == inputNeuronIndex and CurrentSpikeTime > (self.t*second-.1*second):# and CurrentSpikeTime < (self.t*second):
						SpikePreSyn = CurrentSpikeTime
						preSynSpikeFound = True

				# Find SpikePostSyn if exists
				#spikeCollection = M.spikes
				spikeCollection = M.it
				print 'spikeCollection',spikeCollection
				#print 'spikeCollection', spikeCollection
				NumberOfSpikes = len(spikeCollection[0])
				for i in range(NumberOfSpikes):
					#CurrentSpikeNueron = spikeCollection[i][0]
					#CurrentSpikeTime = spikeCollection[i][1]*1000
					CurrentSpikeNueron = spikeCollection[0][i]
					#CurrentSpikeTime = spikeCollection[1][i]*1000
					CurrentSpikeTime = spikeCollection[1][i]*10
					print 'CurrentSpikeNueron',CurrentSpikeNueron,'CurrentSpikeTime',CurrentSpikeTime,'self.t*second',self.t*second

					# exit loop once values below current time elapsed have all been checked
					# Disabled due to spikeCollection not being sorted and causing break too early
					#if CurrentSpikeTime > (self.t*second):
					#	break

					# (self.t*second-.1*second) is to check if in relevant time window below.
					# Note: that may not be a good cut off and I should check it
					# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
					# added (.1+.008)
					if CurrentSpikeNueron == self.neuronIndex and CurrentSpikeTime >= (self.t*second-(.1+.008)*second) and CurrentSpikeTime <= (self.t*second):
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
				print 'neuronIndex', self.neuronIndex, 'inputNeuronIndex',inputNeuronIndex,'self.t', self.t, 'preSynSpikeFound',preSynSpikeFound,'postSynSpikeFound',postSynSpikeFound,'SpikePostSyn', SpikePostSyn, 'SpikePreSyn', SpikePreSyn

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
			print 'self.t',self.t,'SpikePreSyn',SpikePreSyn,'SpikePostSyn',SpikePostSyn,'DeltaSpikeTime',DeltaSpikeTime,'DeltaW',DeltaW
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

		def evaluateClassifier():
			tNorm = self.t - (floor(self.t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
			#totalTimeSteps = 200
			#for t in 1:totalTimeSteps:

			testR = [[1] * numberOfPixels]*dictionaryLongitude
			'''testW = [[0.70106319, 0.00286898, 0.99854988, 0.00329543, 0.99749023, 0.00282395, 0.00237953, 0.29986343, 0.00294557, 0.00333794, 0.99891019, 0.00343499, 0.00283253, 0.70053738, 0.29950164],
			[0.30006903, 0.00455598, 0.99983094, 0.00402913, 0.99792467, 0.00310329, 0.00341758, 0.00273908, 0.70096318, 0.00440027, 0.99921323, 0.00441362, 0.00439755, 0.30008903, 0.70175336],
			[0.00433631, 0.00440846, 0.30120874, 0.00427499, 0.9998187, 0.7008399, 0.00359899, 0.70084357, 0.99970516, 0.0044827, 0.99950441, 0.70234639, 0.00449171, 0.00369084, 0.30203318],
			[0.00364221, 0.00362818, 0.63197624, 0.00505983, 0.9985362, 0.37257988, 0.00489438, 0.99985426, 0.3709176, 0.0038792, 0.99896611, 0.37106075, 0.00546044, 0.00523001, 0.63009559]]'''
			'''testW = [[1, 0, 1,0, 1, 0,0, 0, 0,0, 1, 0,0, 1, 0],
			[0, 0, 1,0, 1, 0,0, 0, 1,0, 1, 0,0, 0, 1],
			[0, 0, 0,0, 1, 1,0, 1, 1,0, 1, 1,0, 0, 0],
			[0, 0, 1,0, 1, 0,0, 1, 0,0, 1, 0,0, 0, 1]]'''
			testW = [[0.7004305, 0.00337241, 0.99925926, 0.00448374, 0.99766777, 0.00371942, 0.00424106, 0.30075355, 0.00287463, 0.00320285, 0.99946073, 0.00240124, 0.00346807, 0.70169285, 0.30123847],
			[0.30042121, 0.00384576, 0.99811889, 0.00313265, 0.9977303, 0.00360853, 0.00338815, 0.0027464, 0.70088053, 0.00438758, 0.99894507, 0.00347196, 0.0026785, 0.30144634, 0.70245531],
			[0.00275161, 0.00460105, 0.300879, 0.00459899, 0.99821505, 0.70091678, 0.00391364, 0.70238888, 0.9995263, 0.00394733, 0.99950849, 0.70174748, 0.00341621, 0.00419915, 0.301429],
			[0.00567632, 0.00332203, 0.63111963, 0.00486383, 0.99876455, 0.37258546, 0.00352507, 0.99790986, 0.37148811, 0.00508246, 0.99973393, 0.37058783, 0.00512713, 0.00412277, 0.63006597]]

			# for each charactor input below test the classfication results
			for neuronIndex in range(dictionaryLongitude):
				# Calculate tauD
				testTauD = self.testTauD
				testTauD = tauDCalc(neuronIndex, testTauD, testW)
				#print 't:',self.t,'neuronIndex',neuronIndex,'main calc tauDCalc()', testTauD

				# Calculate resistance
				testR = self.testR
				testR = resistanceCalc(neuronIndex, testTauD, testR)

				#Sample output for testing
				testR = [[0.010177341248211598, 0.01017372798966064, 0.0065737407129026221, 0.010167162378560104, 0.0045614438315052762, 0.0087232095758120916, 0.010164788112061268, 0.0045655890243921074, 0.0087222022256975579, 0.010164940388151444, 0.0045600840455322393, 0.0087180398636411555, 0.010176198349241663, 0.010165594293192626, 0.0065755845520497269],
				[0.010177341248211598, 0.01017372798966064, 0.0065737407129026221, 0.010167162378560104, 0.0045614438315052762, 0.0087232095758120916, 0.010164788112061268, 0.0045655890243921074, 0.0087222022256975579, 0.010164940388151444, 0.0045600840455322393, 0.0087180398636411555, 0.010176198349241663, 0.010165594293192626, 0.0065755845520497269],
				[0.010177341248211598, 0.01017372798966064, 0.0065737407129026221, 0.010167162378560104, 0.0045614438315052762, 0.0087232095758120916, 0.010164788112061268, 0.0045655890243921074, 0.0087222022256975579, 0.010164940388151444, 0.0045600840455322393, 0.0087180398636411555, 0.010176198349241663, 0.010165594293192626, 0.0065755845520497269],
				[0.010177341248211598, 0.01017372798966064, 0.0065737407129026221, 0.010167162378560104, 0.0045614438315052762, 0.0087232095758120916, 0.010164788112061268, 0.0045655890243921074, 0.0087222022256975579, 0.010164940388151444, 0.0045600840455322393, 0.0087180398636411555, 0.010176198349241663, 0.010165594293192626, 0.0065755845520497269]]

				#print 't:',self.t,'neuronIndex',neuronIndex,'main calc resistanceCalc()', testR

				testId[neuronIndex], self.testIdSpikeFired[neuronIndex], self.testIdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, testId, self.testIdSpikeFired, self.testIdRefractoryPeriod, testR, testW, self.testSpiketimes)
				#print '\ttest\t'
				testUm[neuronIndex] = totalSomaMembranePotential(neuronIndex, testUm, testId, Is, tNorm)
				
				testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex])
				testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex])
				
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
				#print '(tau[neuronIndex][RIndex]*.001)',(tau[neuronIndex][RIndex]*.001)
				#print '((tau[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm)',(((tau[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm)
				#print '(tauM / (tauM - (tau[neuronIndex][RIndex]*.001) ))',(tauM / (tauM - (tau[neuronIndex][RIndex]*.001) ))
				#print 'tauM',tauM
				
				if (tau[neuronIndex][RIndex]*.001 == tauM):
					# avoid a division by 0 issue
					tauM = tauM - .000001

				R[neuronIndex][RIndex] = (((tau[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm) * ((tauM / (tau[neuronIndex][RIndex]*.001) ) ** (tauM / (tauM - (tau[neuronIndex][RIndex]*.001) )))

			return R

		def timePeriodAndRefractoryCalcs():
			neuronIndex = self.neuronIndex
			refractoryPointCounter = self.refractoryPointCounter
			refractoryPointSwitch = self.refractoryPointSwitch
			neuronIndexCounter = self.neuronIndexCounter
			neuronIndexSwitch = self.neuronIndexSwitch
			t = self.t
			
			#t = t + 0.0002
			# changed time step time for testing
			#timeStepInterval = 0.01
			#timeStepInterval = 0.002
			timeStepInterval = 0.001
			#timeStepInterval = 0.0002
			t = t + timeStepInterval
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
				self.lastSpikeInterval = self.t
				WeightChangeCalculation()
				#if self.t == 1.081 or self.t == 2.081 or self.t == 4.081 or self.t == 8.081 or self.t == 12.081:
				'''if self.currentEpoch == 2 or self.currentEpoch == 4 or self.currentEpoch == 8 or self.currentEpoch == 10 or self.currentEpoch == 12:
					print 'self.t:\t',self.t,'\tself.currentEpoch\t',self.currentEpoch,'\tcurrent neuront\t',self.neuronIndex
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
			self.t = t
			self.refractoryPointCounter = refractoryPointCounter
			self.refractoryPointSwitch = refractoryPointSwitch
			self.neuronIndexCounter = neuronIndexCounter			
			self.neuronIndexSwitch = neuronIndexSwitch

		# This network_operation runs the membrane potential calculation function for every milisecond that occurs.
		# The Um output is saved directly to the ADDS V (voltage) records for use with Brian's code.
		@network_operation
		def myoperation():
			print 'neuronIndex', self.neuronIndex, 'self.t', self.t, 'ADDS.t', ADDS.t, 'ADDS.v[0]', ADDS.v
			print 'neuronIndex', self.neuronIndex, 'self.t', self.t, 'ADDS.t', ADDS.t, 'ADDS.IdVZ[0]', ADDS.IdVZ
			#ADDS.IdVZ[0] = 10*mV;ADDS.IdVZ[1] = 10*mV;ADDS.IdVZ[2] = 10*mV;ADDS.IdVZ[3] = 10*mV;

			#print 'neuronIndex', self.neuronIndex, 'self.t', self.t, 'UmM.t',UmM.t
			# Record voltage for whichever neuron is being trained at the time. i.e. neuron 1 for first 300ms.
			# Note: below lline moved to higher code line
			#ADDS[self.neuronIndex][0].V = returnUm(self)
			#ADDS.V[self.neuronIndex] = returnUm(self)
			returnUm(self)
			spikeIntervalCounter = (floor(self.t/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10

			#if refractoryGraphingSwitch == True:
			#ADDS.v2 = 
			#for vCheck in ADDS.v:
			#	if vCheck == -.002*mV:
			#		refractoryGraphingSwitch = True
			for vCheck in range(len(ADDS.v2)):
				if ADDS.v2[vCheck] == 10*mV:
					ADDS.v2[vCheck] = 0*mV

			#ADDS.tNorm[self.neuronIndex] = (self.t - (floor(self.t/spikeIntervalUnformatted) * spikeIntervalUnformatted))*1000*ms
			#print 'clock', myclock.t
			#if ADDS.t[0] >= .101:
			#	print 'bigger'
			#	myclock.reinit()

				#for tIndex in range(len(ADDS.t)):
				#	myclock.reinit()

			'''if ADDS.t[0] >= .0049:
				print 'bigger'
				ADDS.t[0:len(ADDS.t)] = [0*ms]*len(ADDS.t)'''

			## Record spikes
			# removed due to no longer needing
			#if self.UmSpikeFired[self.neuronIndex] == True:
			#	newSpike = [self.neuronIndex, (self.t * ms)]
			#	M.spikes.append(newSpike)

			# classifier performance test
			#evaluateClassifier()

			# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
			# the +1 in (self.totalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
			if spikeIntervalCounter < (self.totalSpikeIntervals+1):
				for neuronIndex in range(dictionaryLongitude):				
					if (self.testUmSpikeFired[neuronIndex] == True) and (spikeIntervalCounter < self.totalSpikeIntervals):
						if (self.correctSpikes[neuronIndex][spikeIntervalCounter] == 1):
							self.truePositiveSpikeResults = self.truePositiveSpikeResults + 1	
							self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = True	
						else:
							self.falsePositiveSpikeResults = self.falsePositiveSpikeResults + 1	
							self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = True

					# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
					# (spikeIntervalCounter-1) is to correct for refractoryPointSwitch occuring after spikeInterval it addresses.
					if self.refractoryPointSwitch == true and (spikeIntervalCounter > 0):
						if self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter-1] == False:
							if (self.correctSpikes[neuronIndex][(spikeIntervalCounter-1)] == 1):
								self.falseNegativeSpikeResults = self.falseNegativeSpikeResults + 1		
							else:
								self.trueNegativeSpikeResults = self.trueNegativeSpikeResults + 1	

						#self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = False

			ADDS.dendriV1[self.neuronIndex] = Id[self.neuronIndex][1]*volt
			ADDS.dendriV2[self.neuronIndex] = Id[self.neuronIndex][2]*volt
			ADDS.dendriV3[self.neuronIndex] = Id[self.neuronIndex][3]*volt
			ADDS.dendriV4[self.neuronIndex] = Id[self.neuronIndex][4]*volt
			ADDS.dendriV5[self.neuronIndex] = Id[self.neuronIndex][5]*volt
			ADDS.dendriV6[self.neuronIndex] = Id[self.neuronIndex][6]*volt
			ADDS.dendriV7[self.neuronIndex] = Id[self.neuronIndex][7]*volt
			ADDS.dendriV8[self.neuronIndex] = Id[self.neuronIndex][8]*volt
			ADDS.dendriV9[self.neuronIndex] = Id[self.neuronIndex][9]*volt
			ADDS.dendriV10[self.neuronIndex] = Id[self.neuronIndex][10]*volt
			ADDS.dendriV11[self.neuronIndex] = Id[self.neuronIndex][11]*volt
			ADDS.dendriV12[self.neuronIndex] = Id[self.neuronIndex][12]*volt
			ADDS.dendriV13[self.neuronIndex] = Id[self.neuronIndex][13]*volt
			ADDS.dendriV14[self.neuronIndex] = Id[self.neuronIndex][14]*volt

			#myclock.tick()
			#defaultclock.tick()

		M = SpikeMonitor(ADDS)
		Mv = StateMonitor(ADDS, 'V', record=True)
		MDendI = StateMonitor(ADDS, 'DendI', record=True)
		MSynI = StateMonitor(ADDS, 'SynI', record=True)
		UmM = StateMonitor(ADDS, 'v', record=True)
		UmM2 = StateMonitor(ADDS, variables=True, record=True)
		UmM3 = StateMonitor(ADDS, 'v2', record=True)
		MdendriV1 = StateMonitor(ADDS, 'dendriV1', record=True)
		MdendriV2 = StateMonitor(ADDS, 'dendriV2', record=True)
		MdendriV3 = StateMonitor(ADDS, 'dendriV3', record=True)
		MdendriV4 = StateMonitor(ADDS, 'dendriV4', record=True)
		MdendriV5 = StateMonitor(ADDS, 'dendriV5', record=True)
		MdendriV6 = StateMonitor(ADDS, 'dendriV6', record=True)
		MdendriV7 = StateMonitor(ADDS, 'dendriV7', record=True)
		MdendriV8 = StateMonitor(ADDS, 'dendriV8', record=True)
		MdendriV9 = StateMonitor(ADDS, 'dendriV9', record=True)
		MdendriV10 = StateMonitor(ADDS, 'dendriV10', record=True)
		MdendriV11 = StateMonitor(ADDS, 'dendriV11', record=True)
		MdendriV12 = StateMonitor(ADDS, 'dendriV12', record=True)
		MdendriV13 = StateMonitor(ADDS, 'dendriV13', record=True)
		MdendriV14 = StateMonitor(ADDS, 'dendriV14', record=True)

		#totalRunTime = 21
		#totalRunTime = 201

		#run(totalRunTime*ms,threads=2, report='text')
		#run(200*ms,threads=2, report='text')
		#run(300*ms,report='text')
		#run(200*ms,report='text')
		run(2000*ms,report='text')
		'''run(10*ms,threads=2, report='text')
		run(10*ms,threads=2, report='text')
		run(10*ms,threads=2, report='text')		'''

		'''store()
		run(100*ms,report='text')
		#myclock.reinit()
		#restore()
		print 'defaultclock.t',defaultclock.t
		restore_initial_state()
		restore()
		ADDS.v[0] = 0*volt;ADDS.v[1] = 0*volt;ADDS.v[2] = 0*volt;ADDS.v[3] = 0*volt
		ADDS.IdVZ[0] = 0*volt;ADDS.IdVZ[1] = 0*volt;ADDS.IdVZ[2] = 0*volt;ADDS.IdVZ[3] = 0*volt
		print 'defaultclock.t',defaultclock.t
		run(100*ms,report='text')
		#myclock.reinit()
		#restore()
		print 'defaultclock.t',defaultclock.t
		restore_initial_state()
		restore()
		ADDS.v[0] = 0*volt;ADDS.v[1] = 0*volt;ADDS.v[2] = 0*volt;ADDS.v[3] = 0*volt
		ADDS.IdVZ[0] = 0*volt;ADDS.IdVZ[1] = 0*volt;ADDS.IdVZ[2] = 0*volt;ADDS.IdVZ[3] = 0*volt		
		print 'defaultclock.t',defaultclock.t
		run(100*ms,report='text')
		#myclock.reinit()'''

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
		
		print 'Final Weights\n0',W[0][:]
		print '1',W[1][:]
		print '2',W[2][:]
		print '3',W[3][:]
		print 'Final Res\n0',R[0][:]
		print '1',R[1][:]
		print '2',R[2][:]
		print '3',R[3][:]		
		#print 'end t',self.t
		print '\n'
		print '+++ Results +++'
		print 'Spike results: TP:\t',self.truePositiveSpikeResults,'\tFP:\t',self.falsePositiveSpikeResults,'\tTN:\t',self.trueNegativeSpikeResults,'\tFN:\t',self.falseNegativeSpikeResults
		print 'totalSpikeIntervalsTested:\t',self.totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
		print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(self.truePositiveSpikeResults, '.1f'))/Decimal(format(self.totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
		print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',(Decimal(format(self.truePositiveSpikeResults, '.1f'))+Decimal(format(self.trueNegativeSpikeResults, '.1f')))/(Decimal(format(self.totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
		print '+++++++++++++++'
		#print 'M.spikes',M.spikes
		#print 'lik spiketimes',spiketimes

		neuronToPlot = 1
		subplot(211)
		#raster_plot(M, title='The gupta network', newfigure=False)
		#raster_plot(M.it, title='The gupta network', newfigure=False)
		plot(M.t/ms, M.i, '.')
		subplot(223)
		'''plot(self.timing, Mv[0] / mV)		
		plot(self.timing, Mv[neuronToPlot] / mV)
		plot(self.timing, Mv[2] / mV)
		plot(self.timing, Mv[3] / mV)'''
		
		#plot(M.t/ms, UmM.v[0])
		#plot(M.t/ms, UmM.v[neuronToPlot])
		#plot(M.t/ms, UmM.v[2])
		#plot(M.t/ms, UmM.v[3])	
		#plot(M.t, UmM.v.T[0])
		#plot(M.t, UmM.v.T[1])
		#plot(M.t, UmM.v.T[2])
		#plot(M.t, UmM.v.T[3])
		#plot(UmM.t, UmM.v.T[0])
		#plot(UmM.t, UmM.v.T[1])
		#plot(UmM.t, UmM.v.T[2])
		#plot(UmM.t, UmM.v.T[3])
		
		#plot(UmM.t, UmM.v.T/mV)
		plot(UmM3.t, UmM3.v2.T/mV)

		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		'''subplot(224)
		plot(MdendriV1.times / ms, MdendriV1[neuronToPlot] / mV)
		plot(MdendriV2.times / ms, MdendriV2[neuronToPlot] / mV)
		plot(MdendriV3.times / ms, MdendriV3[neuronToPlot] / mV)
		plot(MdendriV4.times / ms, MdendriV4[neuronToPlot] / mV)
		plot(MdendriV5.times / ms, MdendriV5[neuronToPlot] / mV)
		plot(MdendriV6.times / ms, MdendriV6[neuronToPlot] / mV)
		plot(MdendriV7.times / ms, MdendriV7[neuronToPlot] / mV)
		plot(MdendriV8.times / ms, MdendriV8[neuronToPlot] / mV)
		plot(MdendriV9.times / ms, MdendriV9[neuronToPlot] / mV)
		plot(MdendriV10.times / ms, MdendriV10[neuronToPlot] / mV)
		plot(MdendriV11.times / ms, MdendriV11[neuronToPlot] / mV)
		plot(MdendriV12.times / ms, MdendriV12[neuronToPlot] / mV)
		plot(MdendriV13.times / ms, MdendriV13[neuronToPlot] / mV)
		plot(MdendriV14.times / ms, MdendriV14[neuronToPlot] / mV)
		xlabel('Time (ms)')
		ylabel('Dendrite Input Current (mV)')'''
		show()

	def __init__(self):
		self.run_model()

def main():
	run_gupta_paper = gupta_paper()

if  __name__ =='__main__':main()
