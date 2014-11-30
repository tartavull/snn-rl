from architecture_further_formulas import *

class gupta_paper:
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	testSpiketimes = spiketimes
	LIK = SpikeGeneratorGroup(15, spiketimes)

	neuronIndex = 0
	neuronCounter = -1
	t = 0.001

	IdSpikeFired = np.array([[False]*numberOfPixels]*dictionaryLongitude); testIdSpikeFired = IdSpikeFired
	IsSpikeFired = np.array([False]*dictionaryLongitude); testIsSpikeFired = IsSpikeFired
	UmSpikeFired = np.array([False]*dictionaryLongitude); testUmSpikeFired = UmSpikeFired

	IdRefractoryPeriod = np.array([[False]*numberOfPixels]*dictionaryLongitude); testIdRefractoryPeriod = IdRefractoryPeriod
	IsRefractoryPeriod = np.array([False]*dictionaryLongitude); testIsRefractoryPeriod = IsRefractoryPeriod
	UmRefractoryPeriod = np.array([False]*dictionaryLongitude); testUmRefractoryPeriod = UmRefractoryPeriod

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

	truePositiveSpikeResults = 0
	falsePositiveSpikeResults = 0
	trueNegativeSpikeResults = 0
	falseNegativeSpikeResults = 0

	totalSpikeIntervals = 12
	correctSpikes = np.array([[1]*totalSpikeIntervals]*dictionaryLongitude);#[[1] * totalSpikeIntervals]*dictionaryLongitude
	correctSpikes[0][3:12] = 0
	correctSpikes[1][0:3] = 0
	correctSpikes[1][6:12] = 0
	correctSpikes[2][0:6] = 0
	correctSpikes[2][9:12] = 0
	correctSpikes[3][0:9] = 0

	testSpikesFiredInInterval = np.array([[False]*(totalSpikeIntervals+1)]*dictionaryLongitude)

	print 'initial Weights\n',W

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
			tNorm = t - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
			spiketimes = self.spiketimes

			# test value Rm
			#Rm = 0.01
			
			# Weight loop
			for WIndex in range(len(W[:][:])):
				if abs(W[neuronIndex][WIndex]) <= 1:
					tauD[neuronIndex][WIndex] = tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)

			#Resistance loop
			for RIndex in range(len(R[:][:])):
				R[neuronIndex][RIndex] = (((tauD[neuronIndex][RIndex]*.001)*neuronFiringThreshold) / Rm) * ((tauM / (tauD[neuronIndex][RIndex]*.001) ) ** (tauM / (tauM - (tauD[neuronIndex][RIndex]*.001) )))

			#Dedritic total post-synaptic current
			#Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(self.neuronIndex)
			Id[neuronIndex], self.IdSpikeFired[neuronIndex], self.IdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, Id, self.IdSpikeFired, self.IdRefractoryPeriod, R, W, spiketimes)

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

			#Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-t/tauS)) + DiracWeightedSum
			Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-tNorm/tauS)) + DiracWeightedSum

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
			#Um[neuronIndex] = -(NeuronInputCoeff - Um2) * (e ** (-t/tauM)) + NeuronInputCoeff
			Um[neuronIndex] = -(NeuronInputCoeff - Um2) * (e ** (-tNorm/tauM)) + NeuronInputCoeff

			# Record membrane potential before refactory period commences
			#ADDS[self.neuronIndex][0].V = Um[neuronIndex] * mV

			#print '\tneuronIndex\t',neuronIndex,'t',t,'SummedDendriteGroup',SummedDendriteGroup,'NeuronInputCoeff',NeuronInputCoeff,'Um[neuronIndex]',Um[neuronIndex]

			# Refractory and spike evaluation: Is
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])
			#self.IsSpikeFired[neuronIndex] = spikeFiredEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])			
			Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(Is[neuronIndex], self.IsSpikeFired[neuronIndex], self.IsRefractoryPeriod[neuronIndex])			

			# Refractory and spike evaluation: Um
			Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])
			#self.UmSpikeFired[neuronIndex] = spikeFiredEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])
			Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(Um[neuronIndex], self.UmSpikeFired[neuronIndex], self.UmRefractoryPeriod[neuronIndex])

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
			timeStepInterval = 0.01
			t = t + timeStepInterval
			self.timing.append(t)
			refractoryPointCounter = refractoryPointCounter + timeStepInterval
			neuronIndexCounter = neuronIndexCounter + timeStepInterval

			# At the end of each spike time interval refractory period is turned off and weight changes
			# are calculated.  Refractory turning off here * is not correct * because it can occur for
			# less than the set refractory period.  I am just using it for the time being for testing.
			refractoryPointSwitch = False
			if refractoryPointCounter >= spikeIntervalUnformatted:
				refractoryPointCounter = 0.001
				refractoryPointSwitch = True
				WeightChangeCalculation()
				#print '*t*',self.t
				#print '***',(self.t == 1.081)
				#print '***2',(self.t == (1081*ms))
				#if self.t == 1.081 or self.t == 2.081 or self.t == 4.081 or self.t == 8.081 or self.t == 12.081:
				'''if self.currentEpoch == 2 or self.currentEpoch == 4 or self.currentEpoch == 8 or self.currentEpoch == 10 or self.currentEpoch == 12:
					print 'self.t:\t',self.t,'\tself.currentEpoch\t',self.currentEpoch,'\tcurrent neuront\t',self.neuronIndex
					print '0',W[0][:]
					print '1',W[1][:]
					print '2',W[2][:]
					print '3',W[3][:]'''

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
			#print 'Um[neuronIndex]',Um[neuronIndex]
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
			#IdResults = []; IdSpikeResults = []; IdRefractResults = [];

			for IdIndex in range(len(IDend[neuronIndex][:])):
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

				Id2 = IDend[neuronIndex][IdIndex]
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
				r = 1
				# for testing .025 below is scaling factor to get spikeCoeff to be at least .011 which is enough to cause spikes.
				SpikeModCoeff = (r*(w*.025))
				#SpikeModCoeff = (r*(w*.012))

				# normalize t to count just time within spike interval.  Formula from article just uses t in a way
				# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
				tNorm = t - (floor(t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
				if Dt <= 0.0:
					DiracFun = 1
					#SpikeModCoeff = (r*w*DiracFun)
					#SpikeModCoeff = .011
					#SpikeModCoeff = .005
					#SpikeModCoeff = 0.25
					SpikeModCoeff = (SpikeModCoeff*DiracFun)
					tauDen = .03	
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff
				else:
					DiracFun = 0
					#SpikeModCoeff = (r*w*DiracFun)
					#SpikeModCoeff = .011
					#SpikeModCoeff = .005
					#SpikeModCoeff = 0.25
					SpikeModCoeff = (SpikeModCoeff*DiracFun)					
					tauDen = .03										
					IDend[neuronIndex][IdIndex] = -(SpikeModCoeff - Id2) * (e ** (-tNorm/tauDen)) + SpikeModCoeff

				# early spike state
				#print 't-.001',(t-.001),'spikeIntervalUnformatted*3',(spikeIntervalUnformatted*3),'((t-.001)==(spikeIntervalUnformatted*3))',((t-.001)==(spikeIntervalUnformatted*3)),'((t-.001)-(spikeIntervalUnformatted*3))',((t-.001)-(spikeIntervalUnformatted*3)),'IdIndex',IdIndex,'neuronIndex',neuronIndex,'Id[neuronIndex][IdIndex] ',Id[neuronIndex][IdIndex],'Dt',Dt,'spikeIntervalUnformatted',spikeIntervalUnformatted,'Dt <= spikeIntervalUnformatted',(Dt <= spikeIntervalUnformatted)

				#print 'before t',t,'IdIndex',IdIndex,'neuronIndex',neuronIndex,'Id[neuronIndex][IdIndex]',Id[neuronIndex][IdIndex],'spikeFired',self.IdSpikeFired[neuronIndex][IdIndex],' ',(self.t - 0.001) % spikeIntervalUnformatted,(self.t - 0.001),spikeIntervalUnformatted
				#print 'before t',t,'IdIndex',IdIndex,'neuronIndex',neuronIndex,'SpikeModCoeff',SpikeModCoeff,'Id[neuronIndex][IdIndex]',Id[neuronIndex][IdIndex],'spikeFired',self.IdSpikeFired[neuronIndex][IdIndex]

				# Refractory and spike evaluation: Id
				#Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex] = refractoryPeriodEvaluation(Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex])
				#self.IdSpikeFired[neuronIndex][IdIndex] = spikeFiredEvaluation(Id[neuronIndex][IdIndex], self.IdSpikeFired[neuronIndex][IdIndex], self.IdRefractoryPeriod[neuronIndex][IdIndex])
				refractResults = refractoryPeriodEvaluation(IDend[neuronIndex][IdIndex], IDendSpikes[neuronIndex][IdIndex], IDendRefract[neuronIndex][IdIndex])
				spikeEvalResults = spikeFiredEvaluation(refractResults[0], refractResults[1], refractResults[2])
				IDend[neuronIndex][IdIndex] = spikeEvalResults[0]
				IDendSpikes[neuronIndex][IdIndex] = spikeEvalResults[1]
				IDendRefract[neuronIndex][IdIndex] = spikeEvalResults[2]
				
				#print 'after t',t,'IdIndex',IdIndex,'neuronIndex',neuronIndex,'SpikeModCoeff',SpikeModCoeff,'Id[neuronIndex][IdIndex]',Id[neuronIndex][IdIndex],'spikeFired',self.IdSpikeFired[neuronIndex][IdIndex],' ',(self.t - 0.001) % spikeIntervalUnformatted,(self.t - 0.001),spikeIntervalUnformatted	
			return [IDend[neuronIndex], IDendSpikes[neuronIndex], IDendRefract[neuronIndex]]

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

		def evaluateClassifier():
			tNorm = self.t - (floor(self.t/spikeIntervalUnformatted) * spikeIntervalUnformatted) 
			#totalTimeSteps = 200
			#for t in 1:totalTimeSteps:

			R = [[1] * numberOfPixels]*dictionaryLongitude
			W = [[0.70106319, 0.00286898, 0.99854988, 0.00329543, 0.99749023, 0.00282395, 0.00237953, 0.29986343, 0.00294557, 0.00333794, 0.99891019, 0.00343499, 0.00283253, 0.70053738, 0.29950164],
			[0.30006903, 0.00455598, 0.99983094, 0.00402913, 0.99792467, 0.00310329, 0.00341758, 0.00273908, 0.70096318, 0.00440027, 0.99921323, 0.00441362, 0.00439755, 0.30008903, 0.70175336],
			[0.00433631, 0.00440846, 0.30120874, 0.00427499, 0.9998187, 0.7008399, 0.00359899, 0.70084357, 0.99970516, 0.0044827, 0.99950441, 0.70234639, 0.00449171, 0.00369084, 0.30203318],
			[0.00364221, 0.00362818, 0.63197624, 0.00505983, 0.9985362, 0.37257988, 0.00489438, 0.99985426, 0.3709176, 0.0038792, 0.99896611, 0.37106075, 0.00546044, 0.00523001, 0.63009559]]

			# for each charactor input below test the classfication results
			for neuronIndex in range(dictionaryLongitude):
				#testId, neuronTestedIdSpikesFired, neuronTestedIdRefractoryPeriod = dentritePostSynapticCurrent()
				testId[neuronIndex], self.testIdSpikeFired[neuronIndex], self.testIdRefractoryPeriod[neuronIndex] = dentritePostSynapticCurrent(neuronIndex, testId, self.testIdSpikeFired, self.testIdRefractoryPeriod, R, W, self.testSpiketimes)
				
				SummedDendriteGroup = sum(testId[neuronIndex])

				NeuronInputCoeff = Rm*(SummedDendriteGroup)

				Um2 = testUm[neuronIndex]
				#NeuronInputCoeff = .011
				#NeuronInputCoeff = .005
				#tauM = .03
				testUm[neuronIndex] = -(NeuronInputCoeff - Um2) * (e ** (-tNorm/tauM)) + NeuronInputCoeff
				#print 'before*1',testUm[neuronIndex]

				testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex] = refractoryPeriodEvaluation(testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex])
				testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex] = spikeFiredEvaluation(testUm[neuronIndex], self.testUmSpikeFired[neuronIndex], self.testUmRefractoryPeriod[neuronIndex])

				#print 'after*1',testUm[neuronIndex]

		# This network_operation runs the membrane potential calculation function for every milisecond that occurs.
		# The Um output is saved directly to the ADDS V (voltage) records for use with Brian's code.
		@network_operation
		def myoperation():
			# Record voltage for whichever neuron is being trained at the time. i.e. neuron 1 for first 300ms.
			# Note: below lline moved to higher code line
			ADDS[self.neuronIndex][0].V = returnUm(self)
			#returnUm(self)
			spikeIntervalCounter = (floor(self.t/spikeIntervalUnformatted) * spikeIntervalUnformatted)*10

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

			# classifier performance test
			evaluateClassifier()

			# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
			# the +1 in (self.totalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
			if spikeIntervalCounter < (self.totalSpikeIntervals+1):
				for neuronIndex in range(dictionaryLongitude):				
					if (self.testUmSpikeFired[neuronIndex] == True) and (spikeIntervalCounter < self.totalSpikeIntervals):
						#print 'neuronIndex\t',neuronIndex,'spikeIntervalCounter\t',spikeIntervalCounter,'\tself.correctSpikes[neuronIndex][spikeIntervalCounter]',self.correctSpikes[neuronIndex][spikeIntervalCounter]
						if (self.correctSpikes[neuronIndex][spikeIntervalCounter] == 1):
							self.truePositiveSpikeResults = self.truePositiveSpikeResults + 1	
							self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = True	
						else:
							self.falsePositiveSpikeResults = self.falsePositiveSpikeResults + 1	
							self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = True

					# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
					# (spikeIntervalCounter-1) is to correct for refractoryPointSwitch occuring after spikeInterval it addresses.
					if self.refractoryPointSwitch == true and (spikeIntervalCounter > 0):
						#print 'rps: neuronIndex\t',neuronIndex,'spikeIntervalCounter\t',spikeIntervalCounter
						if self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter-1] == False:
							if (self.correctSpikes[neuronIndex][(spikeIntervalCounter-1)] == 1):
								self.falseNegativeSpikeResults = self.falseNegativeSpikeResults + 1		
							else:
								self.trueNegativeSpikeResults = self.trueNegativeSpikeResults + 1	

						#self.testSpikesFiredInInterval[neuronIndex][spikeIntervalCounter] = False

					#print 'time\t',self.t, '\t', self.testUmSpikeFired[neuronIndex], '\t',self.correctSpikes[neuronIndex][spikeIntervalCounter],'\t',self.correctSpikeResults

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

		totalRunTime = 21

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
		
		print 'Final Weights\n0',W[0][:]
		print '1',W[1][:]
		print '2',W[2][:]
		print '3',W[3][:]
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
