from architecture_further_formulas import *

class gupta_paper:
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	LIK = SpikeGeneratorGroup(15, spiketimes)

	neuronIndex = 0
	neuronCounter = -1
	t = 0.001

	IdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude)
	IsSpikeFired = np.array([False]*dictionaryLongitude)
	UmSpikeFired = np.array([False]*dictionaryLongitude)

	IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude)
	IsRefractoryPeriod = np.array([False]*dictionaryLongitude)
	UmRefractoryPeriod = np.array([False]*dictionaryLongitude)

	SummedDendriteGroup = 0;
	SynapseToSoma = 0;

	SigmaWDSyn = .0375 * mV
	SigmaWDSynMod = .5

	RmEq = Rm #Rm * mV
	DendR = .0375 * mV
	DendW = 1
	DendDirac = 1	

	currentEpoch = 0

	refractoryPointCounter = 0.0
	refractoryPointSwitch = False
	neuronIndexCounter = 0.0
	neuronIndexSwitch = False

	timing = []

	print 'initial W',W

	def run_model(self):
		
		dictionary = self.dictionary
		spiketimes = self.spiketimes

		eqs = Equations('''
	        V : volt
	        DendI : volt
	        SynI : volt
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

		ADDS = NeuronGroup(N=4, model=eqs,threshold=Vt, reset=Vr)

		def returnUm(self):

			dictionary = self.dictionary
			neuronIndex = self.neuronIndex

			t = self.t
			spiketimes = self.spiketimes
			
			# Weight loop
			for WIndex in range(len(W[:][:])):
				if abs(W[neuronIndex][WIndex]) <= 1:
					tauD[neuronIndex][WIndex] = tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)

			#Resistance loop
			for RIndex in range(len(R[:][:])):
				R[neuronIndex][RIndex] = (((tauD[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm) * ((tauM / (tauD[neuronIndex][RIndex]*.001) ) ** (tauM / (tauM - (tauD[neuronIndex][RIndex]*.001) )))

			#Dedritic total post-synaptic current
			dentritePostSynapticCurrent()

			### Synapse directly to soma ###
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

			Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-t/tauS)) + DiracWeightedSum

			### Soma membrane potential ###
			# Solving for Um in the formula in the article yeilded the below equation
			Um2 = Um[neuronIndex]

			SummedDendriteGroup = sum(Id[neuronIndex])
			SynapseToSoma = Is[neuronIndex]

			self.SummedDendriteGroup = SummedDendriteGroup;
			self.SynapseToSoma = SynapseToSoma;


			#SummedDendriteGroup = 0.00006875
			#SynapseToSoma = 0.00006875

			#NeuronInputCoeff = Rm*(SummedDendriteGroup + SynapseToSoma)
			NeuronInputCoeff = Rm*(SummedDendriteGroup)
			#NeuronInputCoeff = .011
			#tauM = .03
			Um[neuronIndex] = -(NeuronInputCoeff - Um2) * (e ** (-t/tauM)) + NeuronInputCoeff

			#print '\tneuronIndex\t',neuronIndex,'t',t,'SummedDendriteGroup',SummedDendriteGroup,'NeuronInputCoeff',NeuronInputCoeff,'Id:\t',Id

			# Refractory and spike evaluation: Is
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])
			self.IsSpikeFired[neuronIndex] = spikeFiredEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])			

			# Refractory and spike evaluation: Um
			Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])
			self.UmSpikeFired[neuronIndex] = spikeFiredEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])

			#print self.t, 'self.SpikeFired', self.SpikeFired, (spikeIntervalUnformatted), (self.t-0.001) % spikeIntervalUnformatted, ((self.t-0.001) % spikeIntervalUnformatted <= 0.00001), self.refractoryPeriod
			#print 'neuronIndex', neuronIndex, 't', t, 'Um2', Um2, 'SummedDendriteGroup', SummedDendriteGroup, 'SynapseToSoma', SynapseToSoma, 'tauM', tauM
			#print 'neuronIndex', neuronIndex, 't', t, 'weights', W[0][neuronIndex], ' ', W[1][neuronIndex], ' ', W[2][neuronIndex], ' ', W[3][neuronIndex], ' ', W[4][neuronIndex], ' ', W[5][neuronIndex], ' ', W[6][neuronIndex], ' ', W[7][neuronIndex], ' ', W[8][neuronIndex], ' ', W[9][neuronIndex], ' ', W[10][neuronIndex], ' ', W[11][neuronIndex], ' ', W[12][neuronIndex], ' ', W[13][neuronIndex], ' ', W[14][neuronIndex]
			#print 'resistance', R[0][neuronIndex], ' ', R[1][neuronIndex], ' ', R[2][neuronIndex], ' ', R[3][neuronIndex]


			## adjust weights ##
			# weight adjusting needs both presynaptic spike time and post synaptic spike time,
			# therefore it is placed here.
		

			## General STDP learning rule implementation ##
			# To find the SpikePostSyn and SpikePreSyn iterate backward through the spikes
			# list and choose the latest spikes for the relevant neuron.  Those represent
			# the latest pre and post spikes
			# TODO: check this is surely the case




			## include logic for processing if SpikePostSyn or SpikePreSyn were not found.
			## That only occurs in first several seconds though and I don't know it effects
			## the weighting much
			#print 'spikes\t', M.spikes
			#print 'NumberOfSpikes', NumberOfSpikes, '\tSpikePreSyn\t', SpikePreSyn, '\tSpikePostSyn\t',SpikePostSyn

	






			### After all neurons are cycled through the time is iterated by one
			#neuronIndex = neuronIndex + 1
			#if neuronIndex == 4: 
				#neuronIndex = 0
				#t = t + 0.0002
				# changed time step time for testing
				#t = t + 0.002

			refractoryPointCounter = self.refractoryPointCounter
			refractoryPointSwitch = self.refractoryPointSwitch
			neuronIndexCounter = self.neuronIndexCounter
			neuronIndexSwitch = self.neuronIndexSwitch
			
			#t = t + 0.0002
			# changed time step time for testing
			t = t + 0.02
			self.timing.append(t)
			refractoryPointCounter = refractoryPointCounter + 0.02
			neuronIndexCounter = neuronIndexCounter + 0.02

			# At the end of each spike time interval refractory period is turned off and weight changes
			# are calculated.  Refractory turning off here * is not correct * because it can occur for
			# less than the set refractory period.  I am just using it for the time being for testing.
			refractoryPointSwitch = False
			if refractoryPointCounter >= spikeIntervalUnformatted:
				refractoryPointCounter = 0.001
				refractoryPointSwitch = True
				WeightChangeCalculation()

			neuronIndexSwitch = False
			if neuronIndexCounter >= (spikeIntervalUnformatted*3):
				neuronIndexCounter = 0.001
				neuronIndexSwitch = True
				# add epoch every epoch ms duration.  Occurs when a full series of input charactors has
				# occured and a new one starts.  self.neuronIndex == 3 really means it will roll back to
				# self.neuronIndex == 3 after neuronIndex is added
				if self.neuronIndex == 3:
					self.currentEpoch = self.currentEpoch + 1				

			# revision of when new neuron index is used.  new neuron every 300ms
			# .00001 used below due to seeming number slight offset issue perhaps due to number types used
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

			#self.ADDS = ADDS

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

			return SpikeFired

		def dentritePostSynapticCurrent():
			# Solving for Id in the formula in the article yeilded the below equation
			neuronIndex = self.neuronIndex
			t = self.t
			e = math.e

			for IdIndex in range(len(Id[neuronIndex][:])):
				tauDen = tauD[neuronIndex][IdIndex]
				r = R[neuronIndex][IdIndex]
				w = W[neuronIndex][IdIndex]
				# * .001 below is for ms conversion
				#tPreSyn = (spiketimes[IdIndex + (neuronIndex * len(Id)) + (self.currentEpoch*(dictionaryLongitude*numberOfPixels))][1])
				# find most recent spike for neuron

				# set tpresyn initially too far to trigger dirac
				tPreSyn = -t - 1
				for presynInd in range(shape(spiketimes)[0]):
					comparedSpikeTime = spiketimes[presynInd][1]
					'''print 'comparedSpikeTime',comparedSpikeTime
					print 't',t
					print 'spiketimes[presynInd][0]',spiketimes[presynInd][0]
					print 'IdIndex',IdIndex'''
					# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
					# 0.0 to 300 ms, as the input is intended
					# changed <= to < below to see if that works better
					if spiketimes[presynInd][0] == IdIndex and (comparedSpikeTime-spikeIntervalUnformatted) < t:
						## * 0.001 below is conversion to ms
						tPreSyn = comparedSpikeTime
					#print 'IdIndex\t',IdIndex,'\tspiketimes[presynInd][0]\t',spiketimes[presynInd][0],'\tt\t',t,'\tcomparedSpikeTime\t',comparedSpikeTime,'\ttPreSyn\t',tPreSyn
					if (comparedSpikeTime-spikeIntervalUnformatted) > t:
						break

				Id2 = Id[neuronIndex][IdIndex]
				Dt = t - tPreSyn
				# dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
				# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
				# though
				# simplify dirac for testing
				DiracFun = 0

				SpikeModCoeff = (r*w*DiracFun)
				'''print '_______'
				print 'spike time',tPreSyn
				print 't',t
				print 'r',r
				print 'w',w
				print 'tPreSyn',tPreSyn
				print 'DiracFun',DiracFun
				print 'SpikeModCoeff',SpikeModCoeff'''

				# dirac test
				# t in dirac forumula means curent time or last spike time?		
				#if (t > -(Dt/2) and t < (Dt/2)):
				#if Dt <= spikeIntervalUnformatted:
				if Dt <= 0.0:
					DiracFun = 1
					#SpikeModCoeff = (r*w*DiracFun)
					SpikeModCoeff = .011
					#SpikeModCoeff = 0.25
					SpikeModCoeff = (SpikeModCoeff*DiracFun)
					tauDen = .03				
					Id[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-t/tauDen)) + SpikeModCoeff
				else:
					DiracFun = 0
					#SpikeModCoeff = (r*w*DiracFun)
					SpikeModCoeff = .011
					#SpikeModCoeff = 0.25
					SpikeModCoeff = (SpikeModCoeff*DiracFun)					
					tauDen = .03										
					Id[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-t/tauDen)) + SpikeModCoeff

				# early spike state
				#print 't-.001',(t-.001),'spikeIntervalUnformatted*3',(spikeIntervalUnformatted*3),'((t-.001)==(spikeIntervalUnformatted*3))',((t-.001)==(spikeIntervalUnformatted*3)),'((t-.001)-(spikeIntervalUnformatted*3))',((t-.001)-(spikeIntervalUnformatted*3)),'IdIndex',IdIndex,'neuronIndex',neuronIndex,'Id[neuronIndex][IdIndex] ',Id[neuronIndex][IdIndex],'Dt',Dt,'spikeIntervalUnformatted',spikeIntervalUnformatted,'Dt <= spikeIntervalUnformatted',(Dt <= spikeIntervalUnformatted)

				#print 'before t',t,'IdIndex',IdIndex,'neuronIndex',neuronIndex,'Id[neuronIndex][IdIndex]',Id[neuronIndex][IdIndex],'spikeFired',self.IdSpikeFired[neuronIndex][IdIndex],' ',(self.t - 0.001) % spikeIntervalUnformatted,(self.t - 0.001),spikeIntervalUnformatted

				# Refractory and spike evaluation: Id
				Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex] = refractoryPeriodEvaluation(Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex])
				self.IdSpikeFired[neuronIndex][IdIndex] = spikeFiredEvaluation(Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex])
				
				#print 'after t',t,'IdIndex',IdIndex,'neuronIndex',neuronIndex,'Id[neuronIndex][IdIndex]',Id[neuronIndex][IdIndex],'spikeFired',self.IdSpikeFired[neuronIndex][IdIndex],' ',(self.t - 0.001) % spikeIntervalUnformatted,(self.t - 0.001),spikeIntervalUnformatted	

		def WeightChangeCalculation():
			for inputNeuronIndex in range(numberOfPixels):

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
					#print 'checked:\tCurrentSpike[1]',CurrentSpikeTime,'self.t*second',self.t*second,'i',i,'(self.t*second-.1*second)',(self.t*second-.1*second)

					# exit loop once values below current time elapsed have all been checked
					if CurrentSpikeTime > (self.t*second):
						break

					# (self.t*second-.1*second) is to check if in relevant time window below.
					# Note: that may not be a good cut off and I should check it
					if CurrentSpikeNueron == inputNeuronIndex and CurrentSpikeTime > (self.t*second-.1*second):
						SpikePreSyn = CurrentSpikeTime
						preSynSpikeFound = True

					#print 'checked:\tCurrentSpike[1]',CurrentSpikeTime,'SpikePreSyn',SpikePreSyn,'self.t*second',self.t*second,'i',i,'CurrentSpikeNueron',CurrentSpikeNueron,'inputNeuronIndex',inputNeuronIndex,'CurrentSpikeNueron == self.neuronIndex',(CurrentSpikeNueron == self.neuronIndex),'CurrentSpikeTime > (self.t*second-.1*second)',(CurrentSpikeTime > (self.t*second-.1*second)),'(self.t*second-.1*second)',(self.t*second-.1*second)

				# Find SpikePostSyn if exists
				spikeCollection = M.spikes
				NumberOfSpikes = shape(spikeCollection)[0]
				for i in range(NumberOfSpikes):
					CurrentSpikeNueron = spikeCollection[i][0]
					CurrentSpikeTime = spikeCollection[i][1]*1000

					# exit loop once values below current time elapsed have all been checked
					# Disabled due to spikeCollection not being sorted and causing break too early
					#if CurrentSpikeTime > (self.t*second):
					#	break

					# (self.t*second-.1*second) is to check if in relevant time window below.
					# Note: that may not be a good cut off and I should check it
					# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
					if CurrentSpikeNueron == self.neuronIndex and CurrentSpikeTime >= (self.t*second-.1*second) and CurrentSpikeTime <= (self.t*second):
						SpikePostSyn = CurrentSpikeTime	
						postSynSpikeFound = True	

					#print 'checked:\tCurrentSpike[1]',CurrentSpikeTime,'SpikePostSyn',SpikePostSyn,'self.t*second',self.t*second,'i',i,'CurrentSpikeNueron',CurrentSpikeNueron,'inputNeuronIndex',inputNeuronIndex,'CurrentSpikeNueron == self.neuronIndex',(CurrentSpikeNueron == self.neuronIndex),'CurrentSpikeTime > (self.t*second-.1*second)',(CurrentSpikeTime > (self.t*second-.1*second)),'(self.t*second-.1*second)',(self.t*second-.1*second)

				if preSynSpikeFound == False or postSynSpikeFound == False:
					# Todo: watch more for spikes in a relevant time frame but are actually from different synapses and would
					# give different results?  could cause unwanted results? brain would know synapse difference?

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

				# Find DeltaW
				DeltaW = returnDeltaW(SpikePreSyn, SpikePostSyn)  

				# Find new weight
				WOld = W[self.neuronIndex][inputNeuronIndex];

				W[self.neuronIndex][inputNeuronIndex] = returnNewW(WOld, DeltaW);

				#print 'new weight t',self.t,'neuronIndex',self.neuronIndex,'SpikePreSyn',SpikePreSyn,'SpikePostSyn',SpikePostSyn,'inputNeuronIndex',inputNeuronIndex,'WOld',WOld,'w',W[self.neuronIndex][inputNeuronIndex]

		def returnDeltaW(SpikePreSyn, SpikePostSyn):
			DeltaSpikeTime = SpikePreSyn - SpikePostSyn
			# TODO: figure out if DeltaW = 0 is really fine for init value
			DeltaW = 0
			if DeltaSpikeTime < 0:
				DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus*ms)))
				# for testing
				DeltaW = 1
			elif DeltaSpikeTime > 0:
				DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus*ms)))	
				# for testing
				DeltaW = -1

			# testing
			'''print '(t*ms)',(t*ms),'SpikePreSyn',SpikePreSyn,'SpikePostSyn',SpikePostSyn,'DeltaSpikeTime',DeltaSpikeTime,'DeltaW',DeltaW
			if 	(t*ms) - SpikePreSyn > (1*ms):
				DeltaW = .5
			else:
				DeltaW = -.5	'''			

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

		# This network_operation runs the membrane potential calculation function for every milisecond that occurs.
		# The Um output is saved directly to the ADDS V (voltage) records for use with Brian's code.
		@network_operation
		def myoperation():
			# Record voltage for whichever neuron is being trained at the time. i.e. neuron 1 for first 300ms.
			ADDS[self.neuronIndex][0].V = returnUm(self)

			## Record spikes
			'''self.neuronCounter = self.neuronCounter + 1
			if self.neuronCounter == 4:
				self.neuronCounter = 0

			if self.UmSpikeFired[self.neuronCounter] == True:
				newSpike = [self.neuronCounter, (self.t * ms)]
				M.spikes.append(newSpike)'''

			if self.UmSpikeFired[self.neuronIndex] == True:
				newSpike = [self.neuronIndex, (self.t * ms)]
				M.spikes.append(newSpike)

			ADDS[self.neuronIndex][0].dendriV1 = Id[self.neuronIndex][1]
			ADDS[self.neuronIndex][0].dendriV2 = Id[self.neuronIndex][2]
			ADDS[self.neuronIndex][0].dendriV3 = Id[self.neuronIndex][3]
			ADDS[self.neuronIndex][0].dendriV4 = Id[self.neuronIndex][4]
			ADDS[self.neuronIndex][0].dendriV5 = Id[self.neuronIndex][5]
			ADDS[self.neuronIndex][0].dendriV6 = Id[self.neuronIndex][6]
			ADDS[self.neuronIndex][0].dendriV7 = Id[self.neuronIndex][7]
			ADDS[self.neuronIndex][0].dendriV8 = Id[self.neuronIndex][8]
			ADDS[self.neuronIndex][0].dendriV9 = Id[self.neuronIndex][9]
			ADDS[self.neuronIndex][0].dendriV10 = Id[self.neuronIndex][10]
			ADDS[self.neuronIndex][0].dendriV11 = Id[self.neuronIndex][11]
			ADDS[self.neuronIndex][0].dendriV12 = Id[self.neuronIndex][12]
			ADDS[self.neuronIndex][0].dendriV13 = Id[self.neuronIndex][13]
			ADDS[self.neuronIndex][0].dendriV14 = Id[self.neuronIndex][14]

			#print 'neuron\t', self.neuronCounter, '\tADDS[self.neuronCounter][0].Yy\t', ADDS[self.neuronCounter][0].Yy
			#print 'ADDS[self.neuronCounter][0].DendI\t', ADDS[self.neuronCounter][0].DendI

			#print 'ADDS[self.neuronCounter][0].dendriV', ADDS[self.neuronCounter][0].dendriV

		M = SpikeMonitor(ADDS)
		Mv = StateMonitor(ADDS, 'V', record=True)
		MDendI = StateMonitor(ADDS, 'DendI', record=True)
		MSynI = StateMonitor(ADDS, 'SynI', record=True)
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

		totalRunTime = 101

		run(totalRunTime*ms,threads=2, report='text')

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

		#subplot(224)
		#plot(MYy.times / ms, MYy[neuronToPlot] / mV)
		#xlabel('Time (ms)')
		#ylabel('Yy (mV)')
		#print 't',t,'neuronIndex',neuronIndex
		print '0',W[0][:]
		print '1',W[1][:]
		print '2',W[2][:]
		print '3',W[3][:]
		print 'M.spikes',M.spikes
		print 'lik spiketimes',spiketimes

		neuronToPlot = 1
		subplot(211)
		raster_plot(M, title='The gupta network', newfigure=False)
		subplot(223)
		plot(self.timing, Mv[0] / mV)		
		plot(self.timing, Mv[neuronToPlot] / mV)
		plot(self.timing, Mv[2] / mV)
		plot(self.timing, Mv[3] / mV)
		xlabel('Time (ms)')
		ylabel('Membrane Potential (mV)')
		subplot(224)
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
		ylabel('Dendrite Input Current (mV)')
		show()

	def __init__(self):
		self.run_model()

def main():
	run_gupta_paper = gupta_paper()

if  __name__ =='__main__':main()
