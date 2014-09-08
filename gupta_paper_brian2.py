#While reading fninf-08-00006.pdf I copied some useful snippets of code here

#Leaky integrate and fire with refractoriness
G = NeuronGroup ( number_of_neurons,
				'dv/dt = -(v - v_0)/tau_m : volt # membrane potential',
				threshold='v > v_th',
				reset='v = v_0',
				refractory='(t-lastspike) <= 2*ms');

#Random initial values for membrane potential
G.v = 'v_0+randn() *3*mV'

#Spike timming dependant plasticity
S = Synapses ( source_group , target_group,
	'''w: siemens
	dA_source/dt = -A_source/tau_source: siemens (event-driven)
	dA_target/dt = -A_target/tau_target: siemens (event-driven)''',
	pre='''g_post += w
		   A_source += deltaA_source
		   w = clip(w+A_target, 0*siemens, w_max)''',
    post='''A_target += deltaA_target
    		w = clip(w+A_source, 0*siemens, w_max)''')


#Connectivity without self connections
S.connect('i != j')

