from architecture_further_formulas import *

def lateralInhibition(ADDSObj, timeAndRefrac, latInhibSettings):
	# Activate lateral inhibition upon a spike to inhibit other neurons from spiking.  As brian's 
	# Inhibition inhibits before dend input can increase the signal to balance things out.	
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

	for neuronIndex in range(dictionaryLongitude):		
		if neuronIndex != neuronNumberAboutToFire and neuronIndex != neuronNumberFired:
			preliminaryV[neuronIndex] -= abs(latInhibSettings.voltageForInhibition)
			ADDSObj.v[neuronIndex] = preliminaryV[neuronIndex]

	return ADDSObj.v, latInhibSettings