from architecture_further_formulas import *
from weightConnection import *

dictionary = dictionary()
spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
LIK = SpikeGeneratorGroup(15, spiketimes)

'''
class tauDUpdate:
	# Every call to the class iterates the neuronIndex by 1 to return the relevant neuron's values.
	# A W matrix stores the weight values which are updated every time returnTauD is called.
	#Todo: The Resistance loop will be seperate from the Weight loop because it uses results from the Weight loop.
	neuronIndex = 0
	t = 1

	def returnTauD(self):
		neuronIndex = self.neuronIndex
		t = self.t

		# Weight loop
### Updated len(w) below ###
		for WIndex in range(len(W)):
			#print('WeightIndex: ')
			#print(neuronIndex, ' ', WIndex)
			if abs(W[WIndex][neuronIndex]) <= 1:
				tauD[WIndex][neuronIndex] = tauMax - abs(W[WIndex][neuronIndex])*(tauMax-tauMin)

		#Resistance loop
		for RIndex in range(len(R)):
			R[RIndex][neuronIndex] = ((tauD[RIndex][neuronIndex]*neuronFiringThreshold) / Rm) * ((tauM / tauD[RIndex][neuronIndex]) ** (tauM / tauM - tauD[RIndex][neuronIndex]))

		#Dedritic total post-synaptic current
		#%Dpsc = log(-t/stc+((log(-t/stc+tf/stc)*Rd*Wi*heaviside(t-tf))/stc))
		#t = Symbol('t')
		e = math.e
		tauDen = tauD[neuronIndex]
		r = R[neuronIndex]
		w = W[neuronIndex]
		#tPreSyn = tF[NeuronIndex]

		#Below add spike input 1-15
		# in future use numberOfInputNeurons = range(len(Id)) for loop range below?
		### Dentrite ###
		for IdIndex in range(len(Id)):
			tPreSyn = spiketimes[IdIndex + (neuronIndex * len(Id))][0]
			print ('##tPreSyn##')
			print(tPreSyn)
			print(((e ^ (-(-t/tauDen))*r*w*Heaviside(t-tPreSyn))))
			Id[IdIndex][neuronIndex] = (e ^ (-t/tauDen)) + ((e ^ (-(-t/tauDen)+(tPreSyn/tauDen))*r*w*Heaviside(t-tPreSyn))/tauDen)

		### Synapse ###
		for DiracIndex in range(len(Id)):
			tPreSyn = spiketimes[IdIndex + (neuronIndex * len(Id))][0]
			DiracSum = DiracSum + dirac(-t+tPreSyn)
			WSum = WSum + W[WIndex][neuronIndex]
		DiracWeightedSum = (DiracSum + WSum) / len(Id)

		# is the derivative the formula output or Is(t) itself?
		#
		# d(x) at end of formula is just part of integral notation and further derivative does not need
		# to be added in formula below according to http://en.wikipedia.org/wiki/Integral
		Is[neuronIndex] = (e ^ (-t/tauSyn)) + ( (e ^ (-t/tauSyn)) * integral( ( (e ^ (tf[i]/tauSyn)) * (DiracWeightedSum) ) / tauSyn ) )

		### Mem pot ###
		Um[neuronIndex] = (e ^ (-t/tauM)) + ( (e ^ (-t/tauM)) * ( integral( (-(1/(tauM^2) )  ) * (-3 * Rm * tauM * -(e^(tf / tauM)  ) * r * Rm * w * Heaviside(-t+tPreSyn) - Rm * tauM * (integral ( DiracWeightedSum ) ) ) ) ) )

### moved down ###
		### After all neurons are cycled through the time is iterated by one
		neuronIndex = neuronIndex + 1
		if neuronIndex == 4: 
			neuronIndex = 0
			t = t + 1
		self.neuronIndex = neuronIndex
		self.t = t

		####
		## dirac t is:
		## if -Dt/2 < t < Dt/2
		## derivative of a limit as Dt -> 0
		## 1/Dt
		##
		## else
		## limit of 0
		##
		## Seems derivative can be taken the same way as if the Dt -> 0 was not present
		## if wanted I chould check that with someone at some point.
		## since deriv. already includes limit computing something with that should be same
		## references such as http://www.wikihow.com/Do-Derivatives
		####

		#Stub value for now
		return Um[neuronIndex] * mV;
		#return 0 * mV;
		'''

class UmUpdate:
	neuronIndex = 0
	t = 1

	def returnUm(self):
		neuronIndex = self.neuronIndex
		t = self.t

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

			Id[IdIndex][neuronIndex] = ((-Id2**2)/(2*tauDen))+((Id2*r*w)/tPreSyn)

		### Synapse directly to soma ###
		# Solving for Id in the formula in the article yeilded the below equation
		# -Is_1**2/(2*Ts_1) + (Is_1*DiracWeightedSum)/Ts_1
		# That is implemented below
		# To calculate the DiracWeightedSum the spike times with the dirac function applied are multipled by the synapse weight 
		# and summed then divided by the number of synapses for the neuron 
		for DiracIndex in range(len(Is)):
			tPreSyn = spiketimes[DiracIndex + (neuronIndex * len(Is))][0]
			DiracFunctionWithSpikeTimes = (-t+tPreSyn)/tauDen
			DiracWeightedSum = W[WIndex][neuronIndex] * DiracFunctionWithSpikeTimes
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

UmComputation = UmUpdate()

### 100 ms simulation ###
for i in range(4000):
	UmResult = UmComputation.returnUm()
	print 'Time:', UmComputation.t, '\tNeuron Index:', UmComputation.neuronIndex, '\tUmResult:', UmResult

eqs = Equations('''
	  dV/dt  = (-V+ge-gi)/taum : volt
	  dge/dt = -ge/taue        : volt
	  dgi/dt = -gi/taui        : volt
	  ''')
ADDS = NeuronGroup(N=4, model=eqs,threshold=Vt, reset=Vr)

exhitatory2 = weightConnection(LIK, ADDS , 'ge',delay=10*ms,structure='dense')
exhitatory = Connection(LIK, ADDS , 'ge',delay=10*ms,structure='dense')
Wexhitatory = np.random.uniform(10,50,[15,4]) * mV
exhitatory2.connect(LIK,ADDS,Wexhitatory)
Ap = 1 * mV
Am = 1 * mV
stdp=ExponentialSTDP(exhitatory2,taue,taum,Ap,Am,wmax=50 * mV,interactions='all',update='additive')

inhibitory = Connection(ADDS, ADDS , 'gi',delay=5*ms,structure='dense')
#Connect adds layer via lateral inhibitory connections
#the diagonal should be 0 to not auto-inhibate
Winhibitory = np.random.uniform(0,5,[4,4]) * mV
diagonal = np.diag_indices(Winhibitory.shape[0])
Winhibitory[diagonal] = 0;

inhibitory.connect(ADDS,ADDS,Winhibitory)

M = SpikeMonitor(ADDS)
Mv = StateMonitor(ADDS, 'V', record=True)
Mge = StateMonitor(ADDS, 'ge', record=True)
Mgi = StateMonitor(ADDS, 'gi', record=True)

run(10000*ms,threads=2, report='text')

# Present results and logging
if presentResults == True:
	for epochIndex in epochsToPrint:
		reportResultsAndLogging = report_results_and_logging(dictionaryLongitude, epochsToPrint, M, Mv, epochIndex, spikeIntervalUnformatted, dictionary, epochMsDuration)
		reportResultsAndLogging.presenter()	

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
subplot(224)
plot(Mge.times / ms, Mge[neuronToPlot] / mV)
plot(Mgi.times / ms, Mgi[neuronToPlot] / mV)
xlabel('Time (ms)')
ylabel('ge and gi (mV)')
legend(('ge', 'gi'), 'upper right')
show()
