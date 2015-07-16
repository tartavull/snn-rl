import gupta_paper_further_formulas as snnrl
import sys

#sys.argv = ['t1.py', 'test=3', 'reportAccuracy=False']

#snnrl.main()

def snnrlWithParams(*args, **kwargs):
	#sys.argv = ['snnrl.py', args, kwargs]
	#print "\n\nargv:\t",sys.argv[0],sys.argv[1],sys.argv[2],"\n\n"
	return snnrl.main()

#accuracyPerc = snnrlWithParams(["snnrl.py", "evaluateClassifier = True", "accelerateTraining = False", "randomization = 0.75-1.00"])#3#reportAccuracy({'showPlot':'True'})
#accuracyPerc = snnrlWithParams(showPlot=True,randomization='0.75-1.00',reportAccuracy=True)#3#reportAccuracy({'showPlot':'True'})
accuracyPerc = snnrlWithParams()#3#reportAccuracy({'showPlot':'True'})

print accuracyPerc