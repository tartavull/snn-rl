class report_results_and_logging:
	def __init__(self, dictionaryLongitude, epochsToPrint, M, Mv, epochIndex, spikeIntervalUnformatted, dictionary):
		self.epochsToPrint = epochsToPrint
		self.dictionaryLongitude = dictionaryLongitude
		self.M = M
		self.Mv = Mv
		self.epochIndex = epochIndex
		self.spikeIntervalUnformatted = spikeIntervalUnformatted
		self.dictionary = dictionary

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
			SortedMvForEpoch.append(sum(self.Mv[neuronIndex][(self.epochIndex*15):(15+(self.epochIndex*15))]))
		SortedMvForEpoch = sorted(SortedMvForEpoch)
		print ('topVsClosestNtmp: ')
		print SortedMvForEpoch[3]-SortedMvForEpoch[2]
		print ('topVsAllNtmp: ')
		print SortedMvForEpoch[3]-(sum(SortedMvForEpoch[0:(self.dictionaryLongitude-1)])/(self.dictionaryLongitude-1))

	def logger(self, outputFile):
		# The output represents: char presented | epoch number | total neuron membrane voltages in epoch | if neuron had a spike
		if self.epochIndex == 0:
			outputStatement = ['Char\tEpoch\tMembraneVoltage\tSpike\n']
			outputFile.writelines(outputStatement)
		for neuronIndex in range(self.dictionaryLongitude):
			outputStatement = [str(self.dictionary.dictionary[neuronIndex][0]), '\t', str(self.epochIndex), '\t', str(sum(self.Mv[neuronIndex][(self.epochIndex*15):(15+(self.epochIndex*15))])), '\t', str(self.SpikeNumberInEpoch[neuronIndex]), '\n']
			outputFile.writelines(outputStatement)
