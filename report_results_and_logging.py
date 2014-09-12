class report_results_and_logging:
	def __init__(self, dictionaryLongitude, epochsToPrint, M, Mv, epochIndex, spikeIntervalUnformatted, dictionary, epochMsDuration):
		self.epochsToPrint = epochsToPrint
		self.dictionaryLongitude = dictionaryLongitude
		self.M = M
		self.Mv = Mv
		self.epochIndex = epochIndex
		self.spikeIntervalUnformatted = spikeIntervalUnformatted
		self.dictionary = dictionary
		self.epochMsDuration = epochMsDuration

		print 'Epoch number:', self.epochIndex
		SpikeNumberInEpoch = [0] * self.dictionaryLongitude

		for NeuronNumber in range(self.dictionaryLongitude):
			# Spikes per char counter
			for spikeOccurenceTime in self.M[NeuronNumber]:
				if (spikeOccurenceTime >= (self.epochIndex * self.spikeIntervalUnformatted) and spikeOccurenceTime < ((self.epochIndex * self.spikeIntervalUnformatted) + (self.spikeIntervalUnformatted - .001))):
					SpikeNumberInEpoch[NeuronNumber] = SpikeNumberInEpoch[NeuronNumber] + 1

		self.SpikeNumberInEpoch = SpikeNumberInEpoch

	def presenter(self):				
		print ('Avg Percent of Spikes for the Char: ')
		print (sum(self.SpikeNumberInEpoch)/self.dictionaryLongitude)
		print ('Spikes for the Char: ')
		for NeuronNumber in range(self.dictionaryLongitude):
			print self.dictionary.dictionary[NeuronNumber][0], ':', self.SpikeNumberInEpoch[NeuronNumber], ' ',
		print (' ')
		SortedMvForEpoch = []
		for neuronIndex in range(self.dictionaryLongitude):
			SortedMvForEpoch.append(sum(self.Mv[neuronIndex][(self.epochIndex*self.epochMsDuration):(self.epochMsDuration+(self.epochIndex*self.epochMsDuration))]))
		SortedMvForEpoch = sorted(SortedMvForEpoch)
		print ('topVsClosestNtmp: ')
		print SortedMvForEpoch[3]-SortedMvForEpoch[2]
		print ('topVsAllNtmp: ')
		print SortedMvForEpoch[3]-(sum(SortedMvForEpoch[0:(self.dictionaryLongitude-1)])/(self.dictionaryLongitude-1))

	def logger(self, outputFile):
		# The output represents: char presented | epoch number | total positive neuron membrane voltages in epoch | if neuron had a spike
		if self.epochIndex == 0:
			outputStatement = ['Char\tEpoch\tPositiveMembraneVoltage\tSpike\n']
			outputFile.writelines(outputStatement)
		for neuronIndex in range(self.dictionaryLongitude):
			priorMvInEpoch = 0
			totalPositiveMembraneMv = 0
			for MvInEpoch in self.Mv[neuronIndex][(self.epochIndex*self.epochMsDuration):(self.epochMsDuration+(self.epochIndex*self.epochMsDuration))]:
				if MvInEpoch > 0:
					totalPositiveMembraneMv = totalPositiveMembraneMv + MvInEpoch
			outputStatement = [str(self.dictionary.dictionary[neuronIndex][0]), '\t', str(self.epochIndex), '\t', str(totalPositiveMembraneMv), '\t', str(self.SpikeNumberInEpoch[neuronIndex]), '\n']
			outputFile.writelines(outputStatement)
