from architecture_further_formulas import *

def lateralInhibition(ADDSObj, timeAndRefrac, latInhibSettings):
	tNorm = timeAndRefrac.time - (floor((timeAndRefrac.time/.001)*.01) * .1)
	preliminaryV = [None]*dictionaryLongitude
	greatestMembraneVoltage = 0.0
	greatestNeuronNumber = None
	neuronNumberFired = None
	neuronNumberAboutToFire = None
	
	for neuronIndex in range(dictionaryLongitude):
		preliminaryV[neuronIndex] = ADDSObj.prelimV[neuronIndex]
	for neuronIndex in range(dictionaryLongitude):
		for neuronIndex2 in range(len(preliminaryV)):
			if (neuronIndex != neuronIndex2) and (preliminaryV[neuronIndex2] >= greatestMembraneVoltage): 
				#print 'comparison','neuronIndex',neuronIndex,' ',preliminaryV[neuronIndex],'neuronIndex2',neuronIndex2,' ',preliminaryV[neuronIndex2],'greatestMembraneVoltage',greatestMembraneVoltage,preliminaryV[neuronIndex2] >= greatestMembraneVoltage
				greatestMembraneVoltage = preliminaryV[neuronIndex2]
				greatestNeuronNumber = neuronIndex
			if ADDSObj.UmSpikeFired[neuronIndex2] == 1*mV: neuronNumberFired = neuronIndex2

		if (greatestMembraneVoltage >= (10*mV)) and (neuronNumberFired == None):
			latInhibSettings.voltageForInhibition = greatestMembraneVoltage - (10*mV)
			neuronNumberAboutToFire = greatestNeuronNumber
			latInhibSettings.inhibitionReduction = 1.0
		elif (greatestMembraneVoltage < (10*mV)) and (neuronNumberFired == None): 
			latInhibSettings.voltageForInhibition = greatestMembraneVoltage
		else:
			latInhibSettings.voltageForInhibition = 70*mV#(70*mV/latInhibSettings.inhibitionReduction)
			if latInhibSettings.inhibitionReduction < 150: latInhibSettings.inhibitionReduction += 0.01
		#print 'time',timeAndRefrac.time,'neuronIndex',neuronIndex,'neuronIndex2',neuronIndex2,'latInhibSettings.voltageForInhibition',latInhibSettings.voltageForInhibition,(70*mV/latInhibSettings.inhibitionReduction),'latInhibSettings.inhibitionReduction',latInhibSettings.inhibitionReduction

	for neuronIndex in range(dictionaryLongitude):		
		if neuronIndex != neuronNumberAboutToFire and neuronIndex != neuronNumberFired:
			#if latInhibSettings.voltageForInhibition > 0: preliminaryV[neuronIndex] -= latInhibSettings.voltageForInhibition
			#else: preliminaryV[neuronIndex] += latInhibSettings.voltageForInhibition
			preliminaryV[neuronIndex] -= abs(latInhibSettings.voltageForInhibition)
			ADDSObj.v[neuronIndex] = preliminaryV[neuronIndex]

		'''for neuronIndex2 in range(dictionaryLongitude):
			#inhibitionVoltage = ADDSObj.prelimV[neuronIndex2]*latInhibSettings.lateralInhibitScaling
			#if ADDSObj.UmSpikeFired[neuronIndex2] == 1*mV: 
			#	inhibitionVoltage = ADDSObj.inhibitionVoltage[neuronIndex2]#*latInhibSettings.lateralInhibitScaling								
			if neuronIndex != neuronIndex2 and latInhibSettings.voltageForInhibition > 0:
				preliminaryV[neuronIndex] -= inhibitionVoltage'''
		#print 'timeAndRefrac.time',timeAndRefrac.time,'neuronIndex',neuronIndex,'latInhibSettings.voltageForInhibition',latInhibSettings.voltageForInhibition,'ADDSObj.prelimV[neuronIndex]',ADDSObj.prelimV[neuronIndex],'preliminaryV[neuronIndex]',preliminaryV[neuronIndex],'ADDSObj.v[neuronIndex]',ADDSObj.v[neuronIndex],'ADDSObj.UmSpikeFired[neuronIndex]',ADDSObj.UmSpikeFired[neuronIndex]
		#if ADDSObj.UmSpikeFired[neuronIndex] != 1*mV: ADDSObj.v[neuronIndex] = preliminaryV[neuronIndex]

	return ADDSObj.v, latInhibSettings
