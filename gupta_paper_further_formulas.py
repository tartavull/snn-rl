from architecture_further_formulas import *

class gupta_paper:
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	LIK = SpikeGeneratorGroup(15, spiketimes)

	neuronIndex = 0
	neuronCounter = -1
	t = 0.001

	SpikeFired = False

	SummedDendriteGroup = 0;
	SynapseToSoma = 0;

	SigmaWDSyn = .0375 * mV
	SigmaWDSynMod = .5

	RmEq = Rm #Rm * mV
	DendR = .0375 * mV
	DendW = 1
	DendDirac = 1	

	currentEpoch = 0

	refractoryPeriod = False

	def run_model(self):
		
		dictionary = self.dictionary
		spiketimes = self.spiketimes

		# add epoch every epoch ms duration
		if self.t % .1 == 0 & self.neuronIndex == 1:
			self.currentEpoch = self.currentEpoch + 1

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
			for WIndex in range(len(W)):
				if abs(W[WIndex][neuronIndex]) <= 1:
					tauD[WIndex][neuronIndex] = tauMax - abs(W[WIndex][neuronIndex])*(tauMax-tauMin)

			#Resistance loop
			for RIndex in range(len(R)):
				R[RIndex][neuronIndex] = (((tauD[RIndex][neuronIndex]*.001)*neuronFiringThreshold) / Rm) * ((tauM / (tauD[RIndex][neuronIndex]*.001) ) ** (tauM / (tauM - (tauD[RIndex][neuronIndex]*.001) )))

			#Dedritic total post-synaptic current
			# Solving for Id in the formula in the article yeilded the below equation
			e = math.e

			for IdIndex in range(len(Id)):
				tauDen = tauD[IdIndex][neuronIndex]
				r = R[IdIndex][neuronIndex]
				w = W[IdIndex][neuronIndex]
				# * .001 below is for ms conversion
				#tPreSyn = (spiketimes[IdIndex + (neuronIndex * len(Id)) + (self.currentEpoch*(dictionaryLongitude*numberOfPixels))][1])
				# find most recent spike for neuron
				tPreSyn = 0.0
				for presynInd in range(shape(spiketimes)[0]):
					comparedSpikeTime = spiketimes[presynInd][1]*0.001
					if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime < t:
						## * 0.001 below is conversion to ms
						tPreSyn = comparedSpikeTime
					#print 'IdIndex\t',IdIndex,'\tspiketimes[presynInd][0]\t',spiketimes[presynInd][0],'\tt\t',t,'\tcomparedSpikeTime\t',comparedSpikeTime,'\ttPreSyn\t',tPreSyn
					if comparedSpikeTime > t:
						break

				Id2 = Id[IdIndex][neuronIndex]
				Dt = t - tPreSyn
				# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
				# though
				Dt = Dt * 1000
				DiracFun = 1/Dt

				SpikeModCoeff = (r*w*DiracFun)
				#print 'spiketimes',shape(spiketimes)
				#print 'spiketimes',spiketimes[1][:]
				#print spiketimes[1].index(0.1)
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
				if (t > -(Dt/2) and t < (Dt/2)):
					#SpikeModCoeff = (r*w*DiracFun)
					SpikeModCoeff = .011
					SpikeModCoeff = (SpikeModCoeff*DiracFun)
					tauDen = .03					
					'''print 'r',r
					print 'w',w
					print 'tPreSyn',tPreSyn
					print 'DiracFun',DiracFun
					print 'SpikeModCoeff',SpikeModCoeff'''
					Id[IdIndex][neuronIndex] = -(SpikeModCoeff - Id2) * (e ** (-t/tauDen)) + SpikeModCoeff
				else:
					DiracFun = 0
					#SpikeModCoeff = (r*w*DiracFun)
					SpikeModCoeff = .011
					SpikeModCoeff = (SpikeModCoeff*DiracFun)					
					tauDen = .03										
					Id[IdIndex][neuronIndex] = -(SpikeModCoeff - Id2) * (e ** (-t/tauDen)) + SpikeModCoeff

				# Reset threshold.  I assume this current is constricted by this the same as others
				if Id[IdIndex][neuronIndex] >= ActionPotentialThreshold:
					Id[IdIndex][neuronIndex] = Ureset

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
				DiracWeightedSum = DiracWeightedSum + W[WIndex][neuronIndex] * DiracFunctionWithSpikeTimes
			DiracWeightedSum = DiracWeightedSum / len(Id)

			Is2 = Is[neuronIndex]

			Is[neuronIndex] = -(DiracWeightedSum - Is2) * (e ** (-t/tauS)) + DiracWeightedSum

			### Soma membrane potential ###
			# Solving for Um in the formula in the article yeilded the below equation
			Um2 = Um[neuronIndex]

			SummedDendriteGroup = sum(Id[ 0:len(Id) ][neuronIndex])
			SynapseToSoma = Is[neuronIndex]

			self.SummedDendriteGroup = SummedDendriteGroup;
			self.SynapseToSoma = SynapseToSoma;


			#SummedDendriteGroup = 0.00006875
			#SynapseToSoma = 0.00006875

			#NeuronInputCoeff = Rm*(SummedDendriteGroup + SynapseToSoma)
			NeuronInputCoeff = Rm*(SummedDendriteGroup)
			NeuronInputCoeff = .011
			#tauM = .03
			Um[neuronIndex] = -(NeuronInputCoeff - Um2) * (e ** (-t/tauM)) + NeuronInputCoeff


			# Check for spike
			# when identifying a spike, refractoryPeriodEndPoint is used as an alternative condition to refractoryPeriod being
			# over to allow a spike right at that end point

			# Refractory period
			refractoryPeriodEndPoint = (self.t - 0.001) % spikeIntervalUnformatted <= 0.00001			
			if self.SpikeFired == True and (self.t - 0.001) % spikeIntervalUnformatted > 0.00001:
				self.refractoryPeriod = True
				print 'refrac activated!'
			elif self.refractoryPeriod == True and refractoryPeriodEndPoint:
				# 0.001 added above for time step compatibility.   0.00001 instead of 0.0 used for some kind of 
				# modulo computational result offset in python which was producing a value from the modulo calc
				# just slightly over 0
				self.refractoryPeriod = False
				print 'refrac over!'
			#print 'self.refractoryPeriod == True and refractoryPeriodEndPoint', self.refractoryPeriod == True, refractoryPeriodEndPoint, (self.refractoryPeriod == True and refractoryPeriodEndPoint), self.t % (spikeIntervalUnformatted + 0.001) <= 0.00001
			if self.refractoryPeriod == True:
				Um[neuronIndex] = Ureset
				#print self.t

			# Spike firing detection positioned here to allow Ureset value to occur on next timestep above.
			# This allows spike to be recorded in plot instead of value immediately going to Ureset.
			self.SpikeFired = False
			if Um[neuronIndex] >= ActionPotentialThreshold and (self.refractoryPeriod == False or refractoryPeriodEndPoint):
				#Um[neuronIndex] = Ureset
				self.SpikeFired = True				
			#print self.t, 'self.SpikeFired', self.SpikeFired, (spikeIntervalUnformatted), (self.t-0.001) % spikeIntervalUnformatted, ((self.t-0.001) % spikeIntervalUnformatted <= 0.00001), self.refractoryPeriod
			#print 'neuronIndex', neuronIndex, 't', t, 'Um2', Um2#, 'SummedDendriteGroup', SummedDendriteGroup, 'SynapseToSoma', SynapseToSoma, 'tauM', tauM
			#print 'neuronIndex', neuronIndex, 't', t, 'weights', W[0][neuronIndex], ' ', W[1][neuronIndex], ' ', W[2][neuronIndex], ' ', W[3][neuronIndex], ' ', W[4][neuronIndex], ' ', W[5][neuronIndex], ' ', W[6][neuronIndex], ' ', W[7][neuronIndex], ' ', W[8][neuronIndex], ' ', W[9][neuronIndex], ' ', W[10][neuronIndex], ' ', W[11][neuronIndex], ' ', W[12][neuronIndex], ' ', W[13][neuronIndex], ' ', W[14][neuronIndex]
			#print 'resistance', R[0][neuronIndex], ' ', R[1][neuronIndex], ' ', R[2][neuronIndex], ' ', R[3][neuronIndex]


			## adjust weights ##
			# weight adjusting needs both presynaptic spike time and post synaptic spike time,
			# therefore it is placed here.
			# W is weight
			synapseType = 'excitatory'
			# Default 'excitatory' values
			WMin = 0;
			WMax = 1;
			if synapseType == 'inhibitory':
				WMin = -1;
				WMax = 0;				

			## General STDP learning rule implementation ##
			# To find the SpikePostSyn and SpikePreSyn iterate backward through the spikes
			# list and choose the latest spikes for the relevant neuron.  Those represent
			# the latest pre and post spikes
			# TODO: check this is surely the case
			SpikePreSyn = 0.0*usecond
			SpikePostSyn = 0.0*usecond
			SpikeFound = False
			NumberOfSpikes = shape(M.spikes)[0]
			for i in range(NumberOfSpikes):
				CurrentSpike = M.spikes[((NumberOfSpikes - 1) - i)]
				if SpikeFound == True and CurrentSpike[0] == self.neuronIndex:
					# Shift spikes found to pre then post from just pre
					SpikePostSyn = SpikePreSyn
					SpikePreSyn = CurrentSpike[1]
					#exit once both found
					break;
				if SpikeFound == False and CurrentSpike[0] == self.neuronIndex:
					SpikePreSyn = CurrentSpike[1]
					SpikeFound = True
			## include logic for processing if SpikePostSyn or SpikePreSyn were not found.
			## That only occurs in first several seconds though and I don't know it effects
			## the weighting much
			#print 'spikes\t', M.spikes
			#print 'NumberOfSpikes', NumberOfSpikes, '\tSpikePreSyn\t', SpikePreSyn, '\tSpikePostSyn\t',SpikePostSyn

			DeltaSpikeTime = SpikePreSyn - SpikePostSyn
			# TODO: figure out if DeltaW = 0 is really fine for init value
			DeltaW = 0
			if DeltaSpikeTime < 0:
				DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus*ms)))
			elif DeltaSpikeTime > 0:
				DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus*ms)))			

			#print 'DeltaW', DeltaW
			## Relaxation rule implementation ##
			for IdIndex in range(len(Id)):
				WOld = W[IdIndex][neuronIndex];
				if DeltaW >= 0:				
					WNew = WOld + (LearningRate * (DeltaW * (WMax - WOld)))
					
				elif DeltaW < 0:	
					WNew = WOld + (LearningRate * (DeltaW * (WOld - WMin)))
				#print 'WOld', WOld
				#print 'WNew', WNew

				W[IdIndex][neuronIndex] = WNew;
			##

			### After all neurons are cycled through the time is iterated by one
			neuronIndex = neuronIndex + 1
			if neuronIndex == 4: 
				neuronIndex = 0
				#t = t + 0.0002
				# changed time step time for testing
				t = t + 0.002
			self.neuronIndex = neuronIndex
			self.t = t

			#self.ADDS = ADDS

			return Um[neuronIndex] * mV;	

		# This network_operation runs the membrane potential calculation function for every milisecond that occurs.
		# The Um output is saved directly to the ADDS V (voltage) records for use with Brian's code.
		@network_operation
		def myoperation():
			self.neuronCounter = self.neuronCounter + 1
			if self.neuronCounter == 4:
				self.neuronCounter = 0

			ADDS[self.neuronCounter][0].V = returnUm(self)

			## Record spikes
			if self.SpikeFired == True:
				newSpike = [self.neuronCounter, (self.t * ms)]
				M.spikes.append(newSpike)

			ADDS[self.neuronCounter][0].dendriV1 = Id[1][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV2 = Id[2][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV3 = Id[3][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV4 = Id[4][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV5 = Id[5][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV6 = Id[6][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV7 = Id[7][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV8 = Id[8][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV9 = Id[9][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV10 = Id[10][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV11 = Id[11][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV12 = Id[12][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV13 = Id[13][self.neuronCounter]
			ADDS[self.neuronCounter][0].dendriV14 = Id[14][self.neuronCounter]

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

		totalRunTime = 100

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

		neuronToPlot = 1
		subplot(211)
		raster_plot(M, title='The gupta network', newfigure=False)
		subplot(223)
		plot(Mv.times / ms, Mv[0] / mV)		
		plot(Mv.times / ms, Mv[neuronToPlot] / mV)
		plot(Mv.times / ms, Mv[2] / mV)
		plot(Mv.times / ms, Mv[3] / mV)
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
