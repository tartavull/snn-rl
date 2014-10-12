from architecture_further_formulas import *

class gupta_paper:
	dictionary = dictionary()
	spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
	LIK = SpikeGeneratorGroup(15, spiketimes)

	neuronIndex = 0
	neuronCounter = -1
	t = 1
	def run_model(self):
		
		dictionary = self.dictionary
		spiketimes = self.spiketimes

		eqs = Equations('''
		      V : volt
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
				R[RIndex][neuronIndex] = ((tauD[RIndex][neuronIndex]*neuronFiringThreshold) / Rm) * ((tauM / tauD[RIndex][neuronIndex]) ** (tauM / tauM - tauD[RIndex][neuronIndex]))

			#Dedritic total post-synaptic current
			# Solving for Id in the formula in the article yeilded the below equation
			# -Id**2/(2*Td) + Id*Rd*W/Tfd
			# That is implemented below
			e = math.e
			
			for IdIndex in range(len(Id)):
				tauDen = tauD[IdIndex][neuronIndex]
				r = R[IdIndex][neuronIndex]
				w = W[IdIndex][neuronIndex]
				tPreSyn = spiketimes[IdIndex + (neuronIndex * len(Id))][0]
				Id2 = Id[IdIndex][neuronIndex]
				Dt = t - tPreSyn
				DiracFun = 1/Dt

				# dirac test
				# t in dirac forumula means curent time or last spike time?
				if (t > -(Dt/2) and t < (Dt/2)):
					Id[IdIndex][neuronIndex] = ((-Id2**2)/(2*tauDen))+((Id2*r*w)*DiracFun)
				else:
					DiracFun = 0
					Id[IdIndex][neuronIndex] = ((-Id2**2)/(2*tauDen))+((Id2*r*w)*DiracFun)

			### Synapse directly to soma ###
			# Solving for Id in the formula in the article yeilded the below equation
			# -Is_1**2/(2*Ts_1) + (Is_1*DiracWeightedSum)/Ts_1
			# That is implemented below
			# To calculate the DiracWeightedSum the spike times with the dirac function applied are multipled by the synapse weight 
			# and summed then divided by the number of synapses for the neuron 
			DiracWeightedSum = 0
			for DiracIndex in range(len(Is)):
				tPreSyn = spiketimes[DiracIndex + (neuronIndex * len(Is))][0]
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

			Is[neuronIndex] = (-Is2**2)/(2*tauS) + ((Is2*DiracWeightedSum)/tauS)

			### Soma membrane potential ###
			# Solving for Um in the formula in the article yeilded the below equation
			## -Um**2/(2*Tm) + Um*(Rm*SummedDendriteGroupX + Rm*SynapseToSomaX)/Tm
			# That is implemented below
			Um2 = Um[neuronIndex]
			SummedDendriteGroup = sum(Id[ 0:len(Id) ][neuronIndex])
			SynapseToSoma = Is[neuronIndex]

			Um[neuronIndex] = -Um2**2/(2*tauM) + Um2*(Rm*SummedDendriteGroup + Rm*SynapseToSoma)/tauM

			if Um[neuronIndex] == 0: Um[neuronIndex] = Ureset

			### After all neurons are cycled through the time is iterated by one
			neuronIndex = neuronIndex + 1
			if neuronIndex == 4: 
				neuronIndex = 0
				t = t + 1
			self.neuronIndex = neuronIndex
			self.t = t

			return Um[neuronIndex] * mV;	

		# This network_operation runs the membrane potential calculation function for every milisecond that occurs.
		# The Um output is saved directly to the ADDS V (voltage) records for use with Brian's code.
		@network_operation
		def myoperation():
			self.neuronCounter = self.neuronCounter + 1
			if self.neuronCounter == 4:
				self.neuronCounter = 0

			ADDS[self.neuronCounter][0].V = returnUm(self)

		M = SpikeMonitor(ADDS)
		Mv = StateMonitor(ADDS, 'V', record=True)

		totalRunTime = 1000

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

		neuronToPlot = 1
		subplot(211)
		raster_plot(M, title='The gupta network', newfigure=False)
		subplot(223)
		plot(Mv.times / ms, Mv[neuronToPlot] / mV)
		xlabel('Time (ms)')
		ylabel('V (mV)')
		show()

	def __init__(self):
		self.run_model()

def main():
	run_gupta_paper = gupta_paper()

if  __name__ =='__main__':main()
