from architecture_further_formulas import *

def WeightChangeCalculation(neuronIndex, spiketimes, time, negativeWeightReinforcement, positiveWeightReinforcement, M, dendObj):
	for inputNeuronIndex in range(numberOfPixels):
		## General STDP learning rule implementation ##
		# Find SpikePreSyn if it exists
		# Note: could replace this with use of dictionary data structure for lookups if convenient for processing time later
		# Initial values of pre and post with .001 difference used below as estimation of value to use when no post syn spike is 
		# found and weights should decrease.  Based on the formulas the weights would not decrease, as they are intended to, unless
		# a close post syn spike is used in a prior to the pre syn spike position

		SpikePreSyn = time*second
		SpikePostSyn = SpikePreSyn-(.001*second)

		preSynSpikeFound = False
		postSynSpikeFound = False

		spikeCollection = spiketimes
		NumberOfSpikes = shape(spikeCollection)[0]
		for i in range(NumberOfSpikes):
			CurrentSpikeNueron = spikeCollection[i][0]
			CurrentSpikeTime = spikeCollection[i][1]*second

			# exit loop once values below current time elapsed have all been checked
			if CurrentSpikeTime > (time*second):
				break

			# (time*second-.1*second) is to check if in relevant time window below.
			# Note: that may not be a good cut off and I should check it
			if CurrentSpikeNueron == inputNeuronIndex and CurrentSpikeTime > (time*second-.1*second):
				SpikePreSyn = CurrentSpikeTime
				preSynSpikeFound = True

		spikeCollection = M.it
		NumberOfSpikes = len(spikeCollection[0])
		for i in range(NumberOfSpikes):
			CurrentSpikeNueron = spikeCollection[0][i]
			CurrentSpikeTime = spikeCollection[1][i]

			# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
			if CurrentSpikeNueron == neuronIndex and CurrentSpikeTime >= (time*second-(.1+.008)*second) and CurrentSpikeTime <= (time*second):
				SpikePostSyn = CurrentSpikeTime	
				postSynSpikeFound = True	

		if preSynSpikeFound == False or postSynSpikeFound == False:
			SpikePostSyn = SpikePreSyn-(.001*second)
			'''# Todo: watch more for spikes in a relevant time frame but are actually from different synapses and would
			# give different results?  could cause unwanted results? brain would know synapse difference?'''

		# Find DeltaW
		DeltaW, DeltaSpikeTime = returnDeltaW(SpikePreSyn, SpikePostSyn)  

		# Find new weight
		WOld = dendObj[neuronIndex].w[inputNeuronIndex]

		# if statement below skips W change until initial spikes can be formed
		if time>.1: 
			NewW = returnNewW(WOld, DeltaW, DeltaSpikeTime, negativeWeightReinforcement, positiveWeightReinforcement)
			if NewW < 0: NewW = 0
			elif NewW > 1: NewW = 1
			dendObj[neuronIndex].w[inputNeuronIndex] = NewW*volt
		else: 
			dendObj[neuronIndex].w[inputNeuronIndex] = WOld*volt

def returnDeltaW(SpikePreSyn, SpikePostSyn):
	DeltaSpikeTime = SpikePreSyn - SpikePostSyn 
	DeltaSpikeTime = DeltaSpikeTime/second # remove second units
	# In thesis article it explains if DeltaT = 0 then DeltaW = 0
	DeltaW = 0
	# changed below line from that content in the article for a workaround/temp fix
	if DeltaSpikeTime < 0:
		DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus)))
	elif DeltaSpikeTime > 0:
		DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus)))

	return DeltaW, DeltaSpikeTime

def returnNewW(WOld, DeltaW, DeltaSpikeTime, negativeWeightReinforcement, positiveWeightReinforcement):
	WOld = WOld/volt # remove volt unit
	WMin = 0;WMax = 1;
	# Check for inhibition occurence.  DeltaSpikeTime > 0 is used to represent inhibition synapse presence
	if DeltaSpikeTime > 0:
		WMin = -1;WMax = 0;	
	## Relaxation rule implementation ##
	if DeltaW < 0:	
		WNew = WOld + ((LearningRate * (DeltaW * (WOld - WMin)))*negativeWeightReinforcement)
		
	elif DeltaW >= 0:				
		WNew = WOld + ((LearningRate * (DeltaW * (WMax - WOld)))*positiveWeightReinforcement)

	return WNew;
	##	

def tauDCalc(neuronIndex, dendObj):
	# Weight loop
	for WIndex in range(numberOfPixels):
		weight = ((dendObj[neuronIndex].w[WIndex])/volt)
		if abs(weight) <= 1:
			dendObj[neuronIndex].tau[WIndex] = (tauMax - (abs(weight)*(tauMax-tauMin))) * ms# .001 is scaling factor found in dirac method

def resistanceCalc(neuronIndex, dendObj, tauM):
	#Resistance loop
	for RIndex in range(numberOfPixels):				
		if (dendObj[neuronIndex].tau[RIndex]*.001 == tauM*second):
			# avoid a division by 0 issue
			tauM = tauM - .000001

		dendObj[neuronIndex].r[RIndex] = (((((dendObj[neuronIndex].tau[RIndex]*.001)/ms)*neuronFiringThreshold) / Rm) * ((tauM / ((dendObj[neuronIndex].tau[RIndex]/ms)*.001) ) ** (tauM / (tauM - ((dendObj[neuronIndex].tau[RIndex]/ms)*.001) ))))*volt

def diracCalc(IDend, neuronIndex, spiketimes, time, lastSpikeInterval):
	e = math.e
	dend = IDend
	# This method only calculates dirac.  dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
	# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that though

	dendGroup = [None]*len(dend[neuronIndex][:])
	for IdIndex in range(len(dend[neuronIndex][:])):
		# set tpresyn initially too far to trigger dirac
		tPreSyn = time + 1
		for presynInd in range(shape(spiketimes)[0]):
			comparedSpikeTime = spiketimes[presynInd][1]

			if comparedSpikeTime > (lastSpikeInterval + timeAndRefrac.spikeIntervalUnformatted):
				break

			if spiketimes[presynInd][0] == IdIndex and (math.floor(1000*Decimal(format(time, '.3f')))*.001 == math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001):
				tPreSyn = comparedSpikeTime

		Dt = Decimal(format(time, '.8f')) - Decimal(format(tPreSyn, '.8f'))
		
		# simplify dirac until later.  TODO: try more comple dirac   #if (t > -(Dt/2) and t < (Dt/2)):
		# Seems to me what is looked for here is that the post synapse (output neurons) is after the pre synapse (input neurons)
		if Dt >= 0.0:
			DiracFun = 1
			dendGroup[IdIndex] = float(DiracFun)*volt*diracScaling
		else:
			DiracFun = 0
			dendGroup[IdIndex] = float(DiracFun)*volt*diracScaling						

	return dendGroup								