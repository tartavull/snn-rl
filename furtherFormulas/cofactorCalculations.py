from architecture_further_formulas import *

def WeightChangeCalculation(neuronIndex, spiketimes, time, negativeWeightReinforcement, positiveWeightReinforcement, M, dendObj):
	for inputNeuronIndex in range(numberOfPixels):
		## General STDP learning rule implementation ##
		# Find SpikePreSyn if it exists
		# Note: could replace this with use of dictionary data structure for lookups if convenient
		# for processing time later

		# SpikePreSyn not found than make it max distance from SpikePostSyn
		# TODO SpikePreSyn = 0.1 # (Max spike time interval distance)
		# round up to nearest 100ms interval
		SpikePreSyn = time*second#(math.ceil(time*10)*.1)*second
		# SpikePostSyn not found than make it max distance from SpikePreSyn
		# TODO SpikePostSyn = 0.0
		#SpikePostSyn = 0*ms

		# .003 used below as estimation of value to use when no post syn spike is found and weights should
		# decrease.  Based on the formulas the weights would not decrease, as they are intended to, unless
		# a close post syn spike is used in a prior to the pre syn spike position
		SpikePostSyn = SpikePreSyn-(.001*second)#(math.floor(time*10)*.1)*second

		preSynSpikeFound = False
		postSynSpikeFound = False

		spikeCollection = spiketimes
		NumberOfSpikes = shape(spikeCollection)[0]
		for i in range(NumberOfSpikes):
			CurrentSpikeNueron = spikeCollection[i][0]
			CurrentSpikeTime = spikeCollection[i][1]*second

			# exit loop once values below current time elapsed have all been checked
			if CurrentSpikeTime > (time*second):
				#print 'eeww1:\tspiketimes\t',spiketimes,'\tqq\t',CurrentSpikeTime
				break

			# (time*second-.1*second) is to check if in relevant time window below.
			# Note: that may not be a good cut off and I should check it
			# added < time just for testing
			if CurrentSpikeNueron == inputNeuronIndex and CurrentSpikeTime > (time*second-.1*second):# and CurrentSpikeTime < (time*second):
				SpikePreSyn = CurrentSpikeTime
				preSynSpikeFound = True

			#print 'compare~~\t',CurrentSpikeNueron,inputNeuronIndex,CurrentSpikeTime,(time*second-.1*second),preSynSpikeFound

		# Find SpikePostSyn if exists
		#spikeCollection = M.spikes
		spikeCollection = M.it
		#print 'spikeCollection',spikeCollection
		#print 'spikeCollection', spikeCollection
		NumberOfSpikes = len(spikeCollection[0])
		for i in range(NumberOfSpikes):
			#CurrentSpikeNueron = spikeCollection[i][0]
			#CurrentSpikeTime = spikeCollection[i][1]*1000
			CurrentSpikeNueron = spikeCollection[0][i]
			#CurrentSpikeTime = spikeCollection[1][i]*1000
			CurrentSpikeTime = spikeCollection[1][i]#*10
			
			#print 'CurrentSpikeNueron',CurrentSpikeNueron,'neuronIndex',neuronIndex,'CurrentSpikeTime',CurrentSpikeTime,'time*second',time*second

			# exit loop once values below current time elapsed have all been checked
			# Disabled due to spikeCollection not being sorted and causing break too early
			#if CurrentSpikeTime > (time*second):
			#	break

			# (time*second-.1*second) is to check if in relevant time window below.
			# Note: that may not be a good cut off and I should check it
			# * Important difference: CurrentSpikeNueron is compared to self.neuronIndex and not inputNeuronIndex here
			# added (.1+.008)
			if CurrentSpikeNueron == neuronIndex and CurrentSpikeTime >= (time*second-(.1+.008)*second) and CurrentSpikeTime <= (time*second):
				SpikePostSyn = CurrentSpikeTime	
				postSynSpikeFound = True	

		if preSynSpikeFound == False or postSynSpikeFound == False:
			SpikePostSyn = SpikePreSyn-(.001*second)
			'''# Todo: watch more for spikes in a relevant time frame but are actually from different synapses and would
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
				#print 'F F'''

		#SpikePreSyn = 500*ms
		#SpikePostSyn = 700*ms
		
		#print 'neuronIndex', neuronIndex, 'inputNeuronIndex',inputNeuronIndex,'time', time, 'ADDS.v2',ADDS.v2, 'preSynSpikeFound',preSynSpikeFound,'postSynSpikeFound',postSynSpikeFound,'SpikePostSyn', SpikePostSyn, 'SpikePreSyn', SpikePreSyn

		# Find DeltaW
		DeltaW, DeltaSpikeTime = returnDeltaW(SpikePreSyn, SpikePostSyn)  

		# Find new weight
		#WOld = W[neuronIndex][inputNeuronIndex];
		WOld = dendObj[neuronIndex].w[inputNeuronIndex]

		# if statement below skips W change until initial spikes can be formed
		if time>.1: 
			NewW = returnNewW(WOld, DeltaW, DeltaSpikeTime, negativeWeightReinforcement, positiveWeightReinforcement)
			if NewW < 0: NewW = 0
			elif NewW > 1: NewW = 1
			#W[neuronIndex][inputNeuronIndex] = NewW
			dendObj[neuronIndex].w[inputNeuronIndex] = NewW*volt
			#print 'time',time,'NewW > 0',NewW,(NewW > 0),dendObj[neuronIndex].w[inputNeuronIndex]
			#print 'changes::','neuronIndex',neuronIndex,'inputNeuronIndex',inputNeuronIndex,'WOld',WOld,'dendObj[neuronIndex].w[inputNeuronIndex]',dendObj[neuronIndex].w[inputNeuronIndex],'DeltaW',DeltaW,'DeltaSpikeTime',DeltaSpikeTime

			#dendObj[neuronIndex].w[inputNeuronIndex] = W[neuronIndex][inputNeuronIndex]*volt
		else: 
			dendObj[neuronIndex].w[inputNeuronIndex] = WOld*volt
		#print 'DeltaW',DeltaW,'WOld',WOld,'WNew',W[neuronIndex][inputNeuronIndex]
		# Reuse existing values below instead of recalculating
		#somaDirect[neuronIndex].w[inputNeuronIndex] = dendObj[neuronIndex].w[inputNeuronIndex]

def returnDeltaW(SpikePreSyn, SpikePostSyn):
	DeltaSpikeTime = SpikePreSyn - SpikePostSyn 
	DeltaSpikeTime = DeltaSpikeTime/second # remove second units
	# In thesis article it explains if DeltaT = 0 then DeltaW = 0
	DeltaW = 0
	# changed below line from that content in the article for a workaround/temp fix
	#if DeltaSpikeTime < 0:
	if DeltaSpikeTime < 0:
		DeltaW = APlus * (e ** (DeltaSpikeTime / (TauPlus)))
		# for testing
		#DeltaW = 1
	elif DeltaSpikeTime > 0:
		DeltaW = AMinus * (e ** (-DeltaSpikeTime / (TauMinus)))
		# for testing
		#DeltaW = -1

	return DeltaW, DeltaSpikeTime

def returnNewW(WOld, DeltaW, DeltaSpikeTime, negativeWeightReinforcement, positiveWeightReinforcement):
	WOld = WOld/volt # remove volt unit
	WMin = 0;WMax = 1;
	# Check for inhibition occurence.  DeltaSpikeTime > 0 is used to represent inhibition synapse presence
	if DeltaSpikeTime > 0:
		WMin = -1;WMax = 0;	
	#print 'DeltaW', DeltaW
	## Relaxation rule implementation ##
	if DeltaW < 0:	
		WNew = WOld + ((LearningRate * (DeltaW * (WOld - WMin)))*negativeWeightReinforcement)
		
	elif DeltaW >= 0:				
		WNew = WOld + ((LearningRate * (DeltaW * (WMax - WOld)))*positiveWeightReinforcement)

	return WNew;
	##	

def tauDCalc(neuronIndex, dendObj, W):
	# Weight loop
	for WIndex in range(len(W[0][:])):
		if abs(W[neuronIndex][WIndex]) <= 1:
			dendObj[neuronIndex].tau[WIndex] = (tauMax - abs(W[neuronIndex][WIndex])*(tauMax-tauMin)) * ms# .001 is scaling factor found in dirac method

	#return dendObj[neuronIndex].tau

def resistanceCalc(neuronIndex, dendObj, R, tauM):
	#Resistance loop
	for RIndex in range(len(R[0][:])):				
		if (dendObj[neuronIndex].tau[RIndex]*.001 == tauM*second):
			# avoid a division by 0 issue
			tauM = tauM - .000001

		dendObj[neuronIndex].r[RIndex] = (((((dendObj[neuronIndex].tau[RIndex]*.001)/ms)*neuronFiringThreshold) / Rm) * ((tauM / ((dendObj[neuronIndex].tau[RIndex]/ms)*.001) ) ** (tauM / (tauM - ((dendObj[neuronIndex].tau[RIndex]/ms)*.001) ))))*volt

	#return dendObj[neuronIndex].r	

def diracCalc(IDend, neuronIndex, spiketimes, time, lastSpikeInterval):
	e = math.e
	dend = IDend
	# This method only calculates dirac
	# dirac function helps allow or disallow signal to be sent from input neurons through dendrite nodes or not
	# converted Dt to units which seem to make more sense for dirac function used, not sure if right to do that
	# though

	#print ':::range(len(IDend[neuronIndex][:])\t','time\t',time,'\tneuronIndex\t',neuronIndex,'\t',IDend,'\t',IDend[neuronIndex]

	dendGroup = [None]*len(dend[neuronIndex][:])
	for IdIndex in range(len(dend[neuronIndex][:])):
		# set tpresyn initially too far to trigger dirac
		tPreSyn = time + 1
		for presynInd in range(shape(spiketimes)[0]):
			comparedSpikeTime = spiketimes[presynInd][1]

			if comparedSpikeTime > (lastSpikeInterval + timeAndRefrac.spikeIntervalUnformatted):
				break

			#print ':::',comparedSpikeTime,self.lastSpikeInterval,spikeIntervalUnformatted
			##print 'comparedSpikeTime:::\t',comparedSpikeTime,'\t::self.lastSpikeInterval\t',self.lastSpikeInterval,'spikeIntervalUnformatted',spikeIntervalUnformatted,'spiketimes[presynInd][0] == IdIndex',(spiketimes[presynInd][0] == IdIndex),'IdIndex',IdIndex,'comparedSpikeTime > self.lastSpikeInterval',(comparedSpikeTime > self.lastSpikeInterval),'comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted)',(comparedSpikeTime <= (self.lastSpikeInterval + spikeIntervalUnformatted))
			'''# spikeIntervalUnformatted included below to have offset where skiping input starting at .1 becomes 0.0, representing input spikes right from
			# 0.0 to 300 ms, as the input is intended'''
			# checking prior interval for a spike.  This looks for a spike in the prior spike time interval 
			#if spiketimes[presynInd][0] == IdIndex and comparedSpikeTime > lastSpikeInterval and comparedSpikeTime <= (lastSpikeInterval + spikeIntervalUnformatted):
			#print 'IdIndex',IdIndex,'time',time,'comparedSpikeTime',comparedSpikeTime,'time == comparedSpikeTime',(time == comparedSpikeTime),'(Decimal(format(time, .8f)) == Decimal(format(comparedSpikeTime, .8f)))',(Decimal(format(time, '.8f')) == Decimal(format(comparedSpikeTime, '.8f'))),Decimal(format(time, '.8f')),Decimal(format(comparedSpikeTime, '.8f'))
			#if comparedSpikeTime > .89:
			#	print 'compare',(math.floor(1000*Decimal(format(time, '.3f')))*.001 == math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001),math.floor(1000*Decimal(format(time, '.3f')))*.001,math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001
			if spiketimes[presynInd][0] == IdIndex and (math.floor(1000*Decimal(format(time, '.3f')))*.001 == math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001):
				## adding below if statement to avoid constant input stimulus in evaluation, that seem more accurate to how the article did it's evaulation.
				# normalize t to count just time within spike interval.  Formula from article just uses t in a way
				# relative to the spike interval it seems and therefore that approach is applied to the formula here.	
				#tNorm = time - (floor((time/.001)*.01) * .1)
				#print 'tNorm',tNorm
				#if tNorm <= .01:# or evaluationActive == False:
				#if True:
					#tPreSyn = comparedSpikeTime
				tPreSyn = comparedSpikeTime
				#print '++++++++++',time,'+++++++++++','spiketimesnumber',spiketimes[presynInd][0],'neuronIndex',neuronIndex,'spiketimestime',(math.floor(1000*Decimal(format(comparedSpikeTime, '.3f')))*.001),'comparedSpikeTime',comparedSpikeTime

		Dt = Decimal(format(time, '.8f')) - Decimal(format(tPreSyn, '.8f'))
		#print '++++++++Dt++++++++',Dt,'Decimal(format(time, .8f))',Decimal(format(time, '.8f')),'Decimal(format(tPreSyn, .8f))',Decimal(format(tPreSyn, '.8f'))
		# time == lines below are a workaround for initialization values
		if time == 0.121 or time == 0.421 or time == 0.721 or time == 1.021:
			Dt = 1.0
		##print 'neuronIndex',neuronIndex,'tPreSyn000:',tPreSyn,'Dt',Dt
		
		# simplify dirac until later.  TODO: try more comple dirac
		#if (t > -(Dt/2) and t < (Dt/2)):
		# Seems to me what is looked for here is that the post synapse (output neurons) is after the pre synapse (input neurons)
		if Dt >= 0.0:
			if time == 0.121 or time == 0.421 or time == 0.721 or time == 1.021:
				DiracFun = 1.00
			if Dt != 0: 
				# A closer post to pre synapse creates a larger dirac signal for greater weight I'd suppose than
				# a further distance below in DiracFun
				#Dt = Dt*1000
				DiracFun = 1/Dt
			else:
				DiracFun = 1#1000
			DiracFun = 1
			#dendGroup[IdIndex] = float(DiracFun)*volt*.001*self.diracScaling #add .001 as scaling factor.  
			# Scaling factor of .001 could account for 1ms being present instead of .001 as it is computed, difference of
			# 1/1 compared to 1/.001
			dendGroup[IdIndex] = float(DiracFun)*volt*diracScaling
		else:
			DiracFun = 0
			dendGroup[IdIndex] = float(DiracFun)*volt*diracScaling						

	return dendGroup								
