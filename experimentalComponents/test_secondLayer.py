from architecture import *

dictionary = dictionary()
spiketimes = dictionary.spikeTimes(dictionaryLongitude, spikeInterval, spikesPerChar, epochs)
LIK = SpikeGeneratorGroup(15, spiketimes)

eqs = Equations('''
	  dV/dt  = (-V+ge-gi)/taum : volt
	  dge/dt = -ge/taue        : volt
	  dgi/dt = -gi/taui        : volt
	  ''')
ADDS = NeuronGroup(N=4, model=eqs,threshold=Vt, reset=Vr)

exhitatory = Connection(LIK, ADDS , 'ge',delay=10*ms,structure='dense')
Wexhitatory = np.random.uniform(10,50,[15,4]) * mV
exhitatory.connect(LIK,ADDS,Wexhitatory)
Ap = 1 * mV
Am = 1 * mV
stdp=ExponentialSTDP(exhitatory,taue,taum,Ap,Am,wmax=50 * mV,interactions='all',update='additive')


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

run(120 *ms, report='text')

print Wexhitatory

def plotWeight2(data_array, wmax):
	import pylab
	from matplotlib import pyplot, mpl,cm
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt


	num_elem=prod(data_array.shape) #num. of elements to plot
	xpos,ypos=meshgrid(range(data_array.shape[0]),range(data_array.shape[1]))
	xpos=xpos.T.flatten()-0.5 #center bars on integer value of x-axis
	ypos=ypos.T.flatten()-0.5 #center bars on integer value of y-axis
	zpos = zeros(num_elem) #all bars start at z=0
	dx =0.75*ones(num_elem) #width of bars in x-direction
	dy = dx.copy() #width of bars in y-direction (same as x-dir here)
	dz = data_array.flatten() #height of bars from density matrix elements (should use 'real()' if complex)

	nrm=mpl.colors.Normalize(0,max(dz)) #<-- normalize colors to max. data
	#nrm=mpl.colors.Normalize(0,1) #<-- normalize colors to 1
	colors=cm.jet(nrm(dz)) #list of colors for each bar

	#plot figure
	fig = plt.figure(figsize=[6,4])
	ax = Axes3D(fig,azim=-40,elev=70)
	ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors)
	ax.axes.w_xaxis.set_major_locator(IndexLocator(1,-0.5)) #set x-ticks to integers
	ax.axes.w_yaxis.set_major_locator(IndexLocator(1,-0.5)) #set y-ticks to integers
	ax.axes.w_zaxis.set_major_locator(IndexLocator(5,0)) #set z-ticks to integers
	ax.set_zlim3d([0,wmax+1])

	cax,kw=mpl.colorbar.make_axes(ax,shrink=.75,pad=.02) #add colorbar with normalized range
	cb1=mpl.colorbar.ColorbarBase(cax,cmap=cm.jet,norm=nrm)
	plt.show()
	return


def plotWeight(data_array, wmax):
	import pylab
	from matplotlib import pyplot, mpl,cm
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt


	num_elem=prod(data_array.shape) #num. of elements to plot
	xpos,ypos=meshgrid(range(data_array.shape[0]),range(data_array.shape[1]))
	xpos=xpos.T.flatten()-0.5 #center bars on integer value of x-axis
	ypos=ypos.T.flatten()-0.5 #center bars on integer value of y-axis
	zpos = zeros(num_elem) #all bars start at z=0
	dx =0.75*ones(num_elem) #width of bars in x-direction
	dy = dx.copy() #width of bars in y-direction (same as x-dir here)
	dz = data_array.flatten() #height of bars from density matrix elements (should use 'real()' if complex)


	fig = pylab.figure()
	ax = Axes3D(fig)
	ax.axes.w_xaxis.set_major_locator(IndexLocator(1,-0.5)) #set x-ticks to integers
	ax.axes.w_yaxis.set_major_locator(IndexLocator(1,-0.5)) #set y-ticks to integers
	ax.bar(xpos, dz , zs=ypos, zdir='x' , alpha=0.8)

	print xpos,ypos

	plt.show()
	return


plotWeight(Wexhitatory/mV,50) 
