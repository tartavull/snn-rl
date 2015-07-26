from architecture_further_formulas import *
from processHDF5Data import *
from savedWeights import *
import re
import time

def initializeTrainedModelParameters(dendObj):
	# Below values that have been created through prior training are used for w, tau, and r.  Those values are imported
	# through an HD5 file or imported from a raw text version in the savedWeights module

	if loadHD5ForTesting == True:
		hdf5Data = processHDF5Data('simulation.hdf5',permision = "r")
		start = .001 #seconds
		end = 500.0#.5#.13
		weightCollection = np.array([None]*dictionaryLongitude)
		names = [None]*dictionaryLongitude
		for weightIndex in range(dictionaryLongitude): names[weightIndex] = 'weights'+str(weightIndex)
		totalWeightCollection = hdf5Data.incorperateData(weightCollection,names,start,end)
		#print 'len(totalWeightCollection)',len(totalWeightCollection[0]),'shape(totalWeightCollection)',shape(totalWeightCollection[0])
		lastIndex = (len(totalWeightCollection[0])-1)
		testW = []
		for finalWeightIndex in range(dictionaryLongitude):
			testW.extend([totalWeightCollection[finalWeightIndex][lastIndex]])

	else:
		testW = retreiveTestW()

	testW = np.array(testW)
	testTauD = np.array(tauMax - abs(testW) * (tauMax - tauMin))
	tauMScaled = tauM * 1000 # scaling to make the scale match tau
	testR =  np.array(((testTauD * neuronFiringThreshold) / Rm) * ((tauMScaled / testTauD) ** (tauM / (tauMScaled - testTauD))))
	if loadHD5ForTesting == True: del(totalWeightCollection)

	for indexOfDend in range(dictionaryLongitude):
		# TODO see if some of these do not need to be recomputed every time
		dendObj[indexOfDend].w = testW[indexOfDend]*volt*weightScaling # NOTE: weightScaling only included in testing code now, don't be confused with it being in training code!
		dendObj[indexOfDend].tau = testTauD[indexOfDend]*ms # unclear if dividing by tau in ms is e.x. /.02 or /20 but it is assumed to be 20, therefore no ms conversion here
		dendObj[indexOfDend].r = testR[indexOfDend]*mV	# volt unit is cancelled out in the equation anyhow, doesn't matter if it is volt or mV due to being cancelled.  Having mv could cause *.001 that is now wanted
	
	'''print 'Model\'s values\n'
	print 'testW\t',testW
	print 'testTauD\t',testTauD
	print 'testR\r',testR'''

#print 'test start'
#initializeTrainedModelParameters(None)
#print 'test end'

def evaluateClassifierPerf(ADDSObj, testRun):
	# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
	# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
	# the +1 in (totalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
	# * A change was made to the scoring to not count the 0.0-0.1 second period because the input spike generator does not start until .1
	# seconds and the first occurences of output spikes should be monitored looking at .2 seconds to see if any occured in seconds .1-.2 .	
	if timeAndRefrac.refractoryPointSwitch == True:
		if (timeAndRefrac.spikeIntervalCounter >= 2) and (timeAndRefrac.spikeIntervalCounter <= (totalSpikeIntervals+1)):
			for neuronIndex in range(dictionaryLongitude):				
				if (ADDSObj.UmSpikeFired[neuronIndex] == 1*mV):
					if (testRun.correctSpikes[neuronIndex][(timeAndRefrac.spikeIntervalCounter-1)] == 1):
						testRun.truePositiveSpikeResults = testRun.truePositiveSpikeResults + 1	
						print 'TP found\t','self.testSpikeIntervalCounter-1\t',timeAndRefrac.spikeIntervalCounter-1,'neuronIndex\t',neuronIndex
					else:
						testRun.falsePositiveSpikeResults = testRun.falsePositiveSpikeResults + 1	
				elif (ADDSObj.UmSpikeFired[neuronIndex] == 0*mV):
					if (testRun.correctSpikes[neuronIndex][(timeAndRefrac.spikeIntervalCounter-1)] == 1):
						testRun.falseNegativeSpikeResults = testRun.falseNegativeSpikeResults + 1		
					else:
						testRun.trueNegativeSpikeResults = testRun.trueNegativeSpikeResults + 1	
		for neuronIndex in range(dictionaryLongitude):											
			ADDSObj.UmSpikeFired[neuronIndex] = 0*mV

	return ADDSObj.UmSpikeFired, testRun

def evaluateClassifierPerf2(ADDSObj, testRun):
	# Negative results below are only set to be measured after a full spike interval has passed and had the opportunity to have created a spike
	# Only evaluate results for enough epochs to test each char in input (3 spike interv per char * 4 char = 12 spike intervals total)
	# the +1 in (totalSpikeIntervals+1) is to allow a last refractoryPointSwitch triggered negative spike evaluation to occur.
	# * A change was made to the scoring to not count the 0.0-0.1 second period because the input spike generator does not start until .1
	# seconds and the first occurences of output spikes should be monitored looking at .2 seconds to see if any occured in seconds .1-.2 .	

	# Only assign neuron to first window it fires in
	# multiple neurons firing in the same window is not accounted for
	# TODO: properly constrict testing windows for scoring using testingSpikeWindows

	currentSpikeIntervalTime = timeAndRefrac.spikeIntervalCounter
	neuronsWithoutSpikes = 0
	currentWindow = math.floor(Decimal(format((currentSpikeIntervalTime-2), '.1f'))/Decimal(format(testingSpikesPerChar, '.1f')))

	if timeAndRefrac.refractoryPointSwitch == True:
		if (currentSpikeIntervalTime >= 2) and (currentSpikeIntervalTime <= (totalSpikeIntervals+1)):
			for neuronIndex in range(dictionaryLongitude):				
				if (ADDSObj.UmSpikeFired[neuronIndex] == 1*mV):
					if (testRun.firedSpikes[currentSpikeIntervalTime]) == None:
						testRun.firedSpikes[currentSpikeIntervalTime] = np.array([neuronIndex])
					else:
						firedSpikes2 = []
						for elemIndex in range(len(testRun.firedSpikes)): 
							if elemIndex == currentSpikeIntervalTime: 
								newElem = testRun.firedSpikes[elemIndex].tolist()
								newElem.append([neuronIndex])
								firedSpikes2.append(np.array(newElem))
							else: 
								firedSpikes2.append(testRun.firedSpikes[elemIndex])
						testRun.firedSpikes = np.array(firedSpikes2)

					neuronAssignedToWindow = False
					for neuronWindowAssignment in testRun.assignedSpikeWindows:
						if neuronWindowAssignment == currentWindow:
							neuronAssignedToWindow = True

					if neuronAssignedToWindow == False and testRun.assignedSpikeWindows[neuronIndex] == None:
						testRun.assignedSpikeWindows[neuronIndex] = currentWindow

					if testRun.assignedSpikeWindows[neuronIndex] == currentWindow:
						testRun.truePositiveSpikeResults = testRun.truePositiveSpikeResults + 1	
						#print 'TP found\t','self.testSpikeIntervalCounter-1\t',timeAndRefrac.spikeIntervalCounter-1,'neuronIndex\t',neuronIndex
					else:
						testRun.falsePositiveSpikeResults = testRun.falsePositiveSpikeResults + 1
				elif (ADDSObj.UmSpikeFired[neuronIndex] == 0*mV):
					testRun.trueNegativeSpikeResults = testRun.trueNegativeSpikeResults + 1
					neuronsWithoutSpikes += 1

				if neuronsWithoutSpikes == dictionaryLongitude:
					testRun.trueNegativeSpikeResults = testRun.trueNegativeSpikeResults - 1
			correctNeuronFound = False
			#print 'ADDSObj.t',ADDSObj.t,'testRun.firedSpikes', testRun.firedSpikes
			for spikedNeurons in testRun.firedSpikes:
				if len(spikedNeurons) > 1 or spikedNeurons != [None]:
					for spikedNeuronNumber in spikedNeurons:
						if spikedNeuronNumber != None and currentWindow == testRun.assignedSpikeWindows[spikedNeuronNumber]:
							correctNeuronFound = True
			if correctNeuronFound == False:
				testRun.falseNegativeSpikeResults = testRun.falseNegativeSpikeResults + 1
		for neuronIndex in range(dictionaryLongitude):											
			ADDSObj.UmSpikeFired[neuronIndex] = 0*mV

	return ADDSObj.UmSpikeFired, testRun	

def totalCorrectPercentage():
	totalCorrectPercentage = (Decimal(format(testRun.truePositiveSpikeResults, '.1f'))+Decimal(format(testRun.trueNegativeSpikeResults, '.1f')))/(Decimal(format(totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
	return totalCorrectPercentage

def precisionPercentage():
	precisionPercentage = (Decimal(format(testRun.truePositiveSpikeResults, '.1f'))/Decimal(format(totalSpikeIntervals, '.1f')))
	return precisionPercentage	

def OutputEvaluationResults(dendObj, testRun, verbosePrint, evaluateClassifier):
	'''
	If training is performed weights produced are saved to a text file.  The file
	is cleared and formatting rules are applied to put it in a format that can be
	imported into a numpy array.
	Note: this could use code for validating weight files created but that is not writen so far

	Reference: http://stackoverflow.com/questions/22821460/numpy-save-2d-array-to-text-file	
	'''

	if (evaluateClassifier == False):
		weightsFile = open('./furtherFormulas/savedWeights.txt', 'w')
		weightsFile.close()
		weightsFile = open('./furtherFormulas/savedWeights.txt', 'a')
		np.set_printoptions(threshold=np.inf, linewidth=np.inf)  # turn off summarization, line-wrapping
		outputWeightDesc = "["

		for printIndex in range(dictionaryLongitude):
			# Formatting rules applied
			outputWeightDescTemp = ""
			outputWeightDescTemp += np.array2string((dendObj[printIndex].w[:]/volt), separator=', ')
			outputWeightDescTemp = outputWeightDescTemp.replace("0.        ", "0.0")
			outputWeightDescTemp = outputWeightDescTemp.replace("0.]", "0.0]")
			outputWeightDescTemp = outputWeightDescTemp.replace(",  ", ", ")
			outputWeightDescTemp = re.sub(r"0\.\]\,", "0.0]", outputWeightDescTemp)
			outputWeightDescTemp = re.sub(r"[ ]+\]", "]", outputWeightDescTemp)
			outputWeightDescTemp = re.sub(r"[ ]+,", ",", outputWeightDescTemp)
			outputWeightDescTemp = re.sub(r"\[ ", "[", outputWeightDescTemp)
			outputWeightDescTemp = re.sub(r"([0-9]+)\]", r"\1],", outputWeightDescTemp)
			outputWeightDesc += outputWeightDescTemp
			if printIndex < (dictionaryLongitude-1): 
				outputWeightDesc += "\n"
			else:
				outputWeightDesc += "]\n"
				outputWeightDesc = outputWeightDesc.replace("],]", "]]")
				# Remove all brackets with code below because numpy.fromstring used now
				outputWeightDesc = outputWeightDesc.replace("[", "")
				outputWeightDesc = outputWeightDesc.replace("],", "")
				outputWeightDesc = outputWeightDesc.replace("]", "")	
		weightsFile.write(outputWeightDesc)
		weightsFile.close()						

	totalCorrectPercentage = (Decimal(format(testRun.truePositiveSpikeResults, '.1f'))+Decimal(format(testRun.trueNegativeSpikeResults, '.1f')))/(Decimal(format(totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
	if (verbosePrint == True):
		print 'Final Weights\n'
		for printIndex in range(dictionaryLongitude):
			print dendObj[printIndex].w[:]/volt
		print 'Final Tau\n'
		for printIndex in range(dictionaryLongitude):
			print dendObj[printIndex].tau[:]/ms
		print 'Final Res\n'
		for printIndex in range(dictionaryLongitude):
			print dendObj[printIndex].r[:]/mV		
		print '\n'
		print 'testRun.firedSpikes\t',testRun.firedSpikes
		print 'testRun.assignedSpikeWindows\t',testRun.assignedSpikeWindows
		print '\n'			
		print '+++ Results +++'
		print 'Spike results: TP:\t',testRun.truePositiveSpikeResults,'\tFP:\t',testRun.falsePositiveSpikeResults,'\tTN:\t',testRun.trueNegativeSpikeResults,'\tFN:\t',testRun.falseNegativeSpikeResults
		print 'totalSpikeIntervalsTested:\t',totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
		print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(testRun.truePositiveSpikeResults, '.1f'))/Decimal(format(totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
		print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',totalCorrectPercentage
		print '+++++++++++++++'	

	if (evaluateClassifier == True):
		#print ">>>>>>>>>TrueEvalClass<<<<<<<<<<<"
		latestResultsFile = open('./furtherFormulas/latestResults.txt', 'w')
		latestResultsFile.write(str(totalCorrectPercentage))
		latestResultsFile.close()
		time.sleep(20) # Pause to ensure file write is completed before continuing

def printWeights(dendObj):
	print 'Current Weights\n'
	for printIndex in range(dictionaryLongitude):
		print dendObj[printIndex].w[:]/volt	

def writeOptimizationResults(self, testRun, accuracyPerc):
	if (self.optResultsFile != ""):
		optResults = open(self.optResultsFile, 'a')
		optResults.write("\tminWeightsRand\t"+str(self.minWeightsRand)+"\tmaxWeightsRand\t"+str(self.maxWeightsRand)+\
			"\tpositiveWeightReinforcement\t"+str(self.positiveWeightReinforcement)+"\tnegativeWeightReinforcement\t"+str(self.negativeWeightReinforcement))
		optResults.write("\nPrecisionPerc\t"+str(Decimal(format(testRun.truePositiveSpikeResults, '.1f'))/Decimal(format(self.totalSpikeIntervals, '.1f')))+"\t"+str(testRun.truePositiveSpikeResults)+"/"+str(self.totalSpikeIntervals)+"\t(this is the percentage of all true positves that were found)")
		optResults.write('\nAccuracyPerc\t'+str(accuracyPerc)+'\tTP:\t'+str(testRun.truePositiveSpikeResults)+'\tFP:\t'+str(testRun.falsePositiveSpikeResults)+'\tTN:\t'+str(testRun.trueNegativeSpikeResults)+'\tFN:\t'+str(testRun.falseNegativeSpikeResults))
		optResults.close()