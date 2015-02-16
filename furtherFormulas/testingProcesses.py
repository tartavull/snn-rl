from architecture_further_formulas import *

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

	for indexOfDend in range(dictionaryLongitude):
		# TODO see if some of these do not need to be recomputed every time
		dendObj[indexOfDend].w = testW[indexOfDend]*volt*weightScaling
		dendObj[indexOfDend].tau = testTauD[indexOfDend]*ms # unclear if dividing by tau in ms is e.x. /.02 or /20 but it is assumed to be 20, therefore no ms conversion here
		dendObj[indexOfDend].r = testR[indexOfDend]*mV	# volt unit is cancelled out in the equation anyhow, doesn't matter if it is volt or mV due to being cancelled.  Having mv could cause *.001 that is now wanted

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

def OutputEvaluationResults(dendObj, testRun):
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
	print 'Spike results: TP:\t',testRun.truePositiveSpikeResults,'\tFP:\t',testRun.falsePositiveSpikeResults,'\tTN:\t',testRun.trueNegativeSpikeResults,'\tFN:\t',testRun.falseNegativeSpikeResults
	print 'totalSpikeIntervalsTested:\t',totalSpikeIntervals,'\ttotalCharsPresented:\t',dictionaryLongitude
	print 'True positives correct percentage (TP/totalSpikeIntervalsTested):\t',Decimal(format(testRun.truePositiveSpikeResults, '.1f'))/Decimal(format(totalSpikeIntervals, '.1f')),'\t(this is the percentage of all true positves that were found)'
	print 'Total correct percentage (TP+TN/(totalSpikeIntervals*totalCharsPresented)):\t',(Decimal(format(testRun.truePositiveSpikeResults, '.1f'))+Decimal(format(testRun.trueNegativeSpikeResults, '.1f')))/(Decimal(format(totalSpikeIntervals, '.1f'))*Decimal(format(dictionaryLongitude, '.1f')))
	print '+++++++++++++++'	