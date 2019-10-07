def printOutputForTesting(neuronIndex, timeAndRefrac, dendObj):
	print( 'neuronIndex',neuronIndex,'time',timeAndRefrac.time)
	print( 'dendObj.dirac',dendObj[0].dirac,dendObj[1].dirac,dendObj[2].dirac,dendObj[3].dirac)
	print( 'dendObj.v',dendObj[0].v,dendObj[1].v,dendObj[2].v,dendObj[3].v)
	print( 'sum(dendObj[neuronIndex].v[:])',sum(dendObj[neuronIndex].v[:]))
	print( 'Prelim Weights\n0',dendObj[0].w[:])
	print( '1',dendObj[1].w[:])
	print( '2',dendObj[2].w[:])
	print( '3',dendObj[3].w[:])
	print( 'Prelim Tau\n0',dendObj[0].tau[:])
	print( '1',dendObj[1].tau[:])
	print( '2',dendObj[2].tau[:])
	print( '3',dendObj[3].tau[:]					)
	print( 'Prelim Res\n0',dendObj[0].r[:])
	print( '1',dendObj[1].r[:])
	print( '2',dendObj[2].r[:])
	print( '3',dendObj[3].r[:])

def sortVersionOutForTesting(neuronIndex, timeAndRefrac, ADDSObj, dendObj):
	print( 'ADDSObj.v','time',timeAndRefrac.time,ADDSObj.v[0],ADDSObj.v[1],ADDSObj.v[2],ADDSObj.v[3])
	print( 'sum(ADDSObj.DendI)',(ADDSObj.DendI[0]),(ADDSObj.DendI[1]),(ADDSObj.DendI[2]),(ADDSObj.DendI[3]))
	print( 'sum(ADDSObj.SynI)',sum(ADDSObj.SynI[0]),sum(ADDSObj.SynI[1]),sum(ADDSObj.SynI[2]),sum(ADDSObj.SynI[3])						)

def printEvalOutputForTesting(neuronIndex, timeAndRefrac, ADDSObj, dendObj):
	'''print 'dendObj.dirac','time',timeAndRefrac.time,'\t',dendObj[0].dirac,dendObj[1].dirac,dendObj[2].dirac,dendObj[3].dirac
	print 'dendObj.v',dendObj[0].v,dendObj[1].v,dendObj[2].v,dendObj[3].v
	print 'sum(dendObj[neuronIndex])',sum(dendObj[0].v[:]),sum(dendObj[1].v[:]),sum(dendObj[2].v[:]),sum(dendObj[3].v[:])							
	print 'testSomaDirect.v',testSomaDirect.v[0],testSomaDirect.v[1],testSomaDirect.v[2],testSomaDirect.v[3]
	print 'sum(testSomaDirect[neuronIndex])',sum(testSomaDirect.v[0]),sum(testSomaDirect.v[1]),sum(testSomaDirect.v[2]),sum(testSomaDirect.v[3])							
	print 'testSomaDirect.summedWandDirac',testSomaDirect.summedWandDirac[0],testSomaDirect.summedWandDirac[1],testSomaDirect.summedWandDirac[2],testSomaDirect.summedWandDirac[3]					'''
	print( 'ADDSObj.v',ADDSObj.v[0],ADDSObj.v[1],ADDSObj.v[2],ADDSObj.v[3])
	print( 'dendObj.dirac','time',timeAndRefrac.time,'\t',dendObj[0].dirac,dendObj[1].dirac,dendObj[2].dirac,dendObj[3].dirac)
	print( 'dendObj.v',dendObj[0].v,dendObj[1].v,dendObj[2].v,dendObj[3].v)
	print( 'v',dendObj[0].v,'\t',dendObj[1].v)
	print( 'r',dendObj[0].r,'\t',dendObj[1].r)
	print( 'w',' time',timeAndRefrac.time,'\t',dendObj[0].w,'\t',dendObj[1].w,'\t',dendObj[2].w,'\t',dendObj[3].w)
	print( 'dirac',dendObj[0].dirac,'t',dendObj[1].dirac)
	print( 'tau',dendObj[0].tau,'t',dendObj[1].tau)
	print( 'units',(-dendObj[1].v[0]/mV),((dendObj[1].r[0])/volt),((dendObj[1].w[0])/volt),(dendObj[1].dirac[0]/volt),(dendObj[1].tau[0]/second))
	dendObj[1].vTest[0] = (((-dendObj[1].v[0]/mV)+(((dendObj[1].r[0])/volt)*((dendObj[1].w[0])/volt)*(dendObj[1].dirac[0]/volt)))/(dendObj[1].tau[0]/second))*volt#((-dendObj[1].v[0]/volt)+((dendObj[1].r[0])/volt))/(dendObj[1].tau[0]/ms)*volt
	dendObj[1].vTest2[0] = -dendObj[1].v[0]
	dendObj[1].vTest3[0] = ((dendObj[1].r[0]/volt)*(dendObj[1].w[0]/volt)*(dendObj[1].dirac[0]))/(dendObj[1].tau[0]/second)
	print( 'vTest',dendObj[1].vTest[0])
	print( 'vTest2',dendObj[2].vTest2[0],dendObj[1].v[0],dendObj[1].v[1],dendObj[1].v[2],dendObj[1].v[3],dendObj[1].v[4])
	print( 'vTest3',dendObj[3].vTest3[0])
	#print 'rs equ',(-dendObj[0].v+((dendObj[0].r)*(dendObj[0].w)*dendObj[0].dirac)/(dendObj[0].tau)),'\t',(-dendObj[1].v+((dendObj[1].r)*(dendObj[1].w)*dendObj[1].dirac)/(dendObj[1].tau))
	print( 'sum(ADDSObj.DendI)',(ADDSObj.DendI[0]),(ADDSObj.DendI[1]),(ADDSObj.DendI[2]),(ADDSObj.DendI[3]))
	print( 'sum(ADDSObj.SynI)',sum(ADDSObj.SynI[0]),sum(ADDSObj.SynI[1]),sum(ADDSObj.SynI[2]),sum(ADDSObj.SynI[3]))
	print( 'ADDSObj.UmSpikeFired',ADDSObj.UmSpikeFired[0],ADDSObj.UmSpikeFired[1],ADDSObj.UmSpikeFired[2],ADDSObj.UmSpikeFired[3])
	#print 'Prelim Weights\n0',dendObj[0].w[:]
	#print '1',dendObj[1].w[:]
	#print '2',dendObj[2].w[:]
	#print '3',dendObj[3].w[:]
	#print 'Prelim Tau\n0',dendObj[0].tau[:]
	#print '1',dendObj[1].tau[:]
	#print '2',dendObj[2].tau[:]
	#print '3',dendObj[3].tau[:]
	#print 'Prelim Res\n0',dendObj[0].r[:]
	#print '1',dendObj[1].r[:]
	#print '2',dendObj[2].r[:]
	#print '3',dendObj[3].r[:]
	print('____________________')
