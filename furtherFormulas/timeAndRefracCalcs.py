from architecture_further_formulas import *

def timePeriodAndRefractoryCalcs(timeAndRefrac):
	timeAndRefrac.time = timeAndRefrac.time + timeAndRefrac.timeStepInterval
	timeAndRefrac.refractoryPointCounter = timeAndRefrac.refractoryPointCounter + timeAndRefrac.timeStepInterval	
	tT = (Decimal(format(timeAndRefrac.time, '.1f')))
	sIU = (Decimal(format(timeAndRefrac.spikeIntervalUnformatted, '.1f')))
	timeAndRefrac.spikeIntervalCounter = int(tT/sIU)

	# if statement below corrects for offset of spiketimes starting at .1 sec.
	if timeAndRefrac.time >= .1:
		timeAndRefrac.neuronIndexCounter = timeAndRefrac.neuronIndexCounter + timeAndRefrac.timeStepInterval

	timeAndRefrac.refractoryPointSwitch = False
	if timeAndRefrac.refractoryPointCounter >= timeAndRefrac.spikeIntervalUnformatted:
		timeAndRefrac.refractoryPointCounter = 0.000
		timeAndRefrac.refractoryPointSwitch = True
		timeAndRefrac.lastSpikeInterval = timeAndRefrac.time-timeAndRefrac.timeStepInterval

	return timeAndRefrac

def checkForResets(neuronIndex, ADDSObj, dendObj, somaDirectObj, timeAndRefrac):
	if ADDSObj.v2[neuronIndex] == 10*mV:
		ADDSObj.v2[neuronIndex] = 0*mV						

	elif ADDSObj.v[neuronIndex] != -0.002 * mV:
		ADDSObj.v2[neuronIndex] = ADDSObj.v[neuronIndex]

	if somaDirectObj.v2[neuronIndex] == 10*mV:
		somaDirectObj.v2[neuronIndex] = 0*mV	

	if (timeAndRefrac.beginRefrac[neuronIndex] == 1*mV) and (timeAndRefrac.relRefracTimeActivated[neuronIndex] == False):
		timeAndRefrac.beginRefrac[neuronIndex] = 0*mV
		timeAndRefrac.absRefracTime[neuronIndex] = timeAndRefrac.time
		timeAndRefrac.relRefracTime[neuronIndex] = timeAndRefrac.time
		timeAndRefrac.absRefracTimeActivated[neuronIndex] = True
		timeAndRefrac.relRefracTimeActivated[neuronIndex] = True
	# Absolute refrac
	if timeAndRefrac.absRefracTimeActivated[neuronIndex] == True:
		if timeAndRefrac.time >= (timeAndRefrac.absRefracTime[neuronIndex] + timeAndRefrac.absRefracTimeDuration):
			timeAndRefrac.absRefracTimeActivated[neuronIndex] = False
	# Relative refrac
	if timeAndRefrac.relRefracTimeActivated[neuronIndex] == True:
		if timeAndRefrac.time >= (timeAndRefrac.relRefracTime[neuronIndex] + timeAndRefrac.relRefracTimeDuration):
			timeAndRefrac.relRefracTimeActivated[neuronIndex] = False

	return ADDSObj.v2, somaDirectObj.v2, timeAndRefrac