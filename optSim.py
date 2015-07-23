from hyperopt import fmin, tpe, hp
#import gupta_paper_further_formulas as snnrl
import sys
from subprocess import * 
from decimal import Decimal

simFile = "gupta_paper_further_formulas.py"
optResultsFile = "optimizationResults.txt"
trainingRunTime = 33000#1100#1100
testingRunTime = 7900
minRandVals = ['minRand', 0.5, 0.75]
maxRandVals = ['maxRand', 0.76, 1.0]
posReinfVals = ['posReinfVal', 0.0, 5.0]
negReinfVals = ['negReinfVal', 0.0, 5.0]
space = [hp.uniform(minRandVals[0],minRandVals[1],minRandVals[2]), hp.uniform(maxRandVals[0],maxRandVals[1],maxRandVals[2]), \
         hp.uniform(posReinfVals[0],posReinfVals[1],posReinfVals[2]), hp.uniform(negReinfVals[0],negReinfVals[1],negReinfVals[2])]
evaluationRuns = 6;

# Clear and open output file
optResults = open(optResultsFile, 'w')
optResults.write("Parameters used in run:\n")
optResults.write("trainingRunTime:\t"+str(trainingRunTime)+"\n")
optResults.write("testingRunTime:\t\t"+str(testingRunTime)+"\tNote: testingRunTime parameter passing not currently working\n")
optResults.write(str(minRandVals)+"\n"+str(maxRandVals)+"\n"+str(posReinfVals)+"\n"+str(negReinfVals)+"\n")
optResults.close()

def snnrlWithParams(kwargs):
	#sys.argv = ['snnrl.py', args, kwargs]
    kwargs = kwargs.split(", ")
    print "out:\n\n",kwargs
    simRun = Popen(["python", simFile]+kwargs, stdout=PIPE)
    simRun.wait()
    results = simRun.communicate()[0].strip()
    print "\n\nresults:\n"+results+"\n\n"
    if str(results) != '': results = Decimal(results, '.1f')
    return results

def objective(args):
    minRand, maxRand, posReinfVal, negReinfVal = args
    # train
    randVal = str(minRand)+'-'+str(maxRand)
    snnrlWithParams("evaluateClassifier=False, runTime="+str(trainingRunTime)+", standardPrint=False, verbosePrint=False, randomization="+str(randVal)+", posReinf="+str(posReinfVal)+", negReinf="+str(negReinfVal))
    # test
    results = snnrlWithParams("evaluateClassifier=True, optResultsFile=\'"+str(optResultsFile)+"\', standardPrint=False, verbosePrint=False, randomization="+str(randVal)+", posReinf="+str(posReinfVal)+", negReinf="+str(negReinfVal))
    accuracyPerc = Decimal(1, '.1f') - results
    print "\naccPer:\t",accuracyPerc
    return {'loss': accuracyPerc}

best = fmin(objective,
    space,
    algo=tpe.suggest,
    max_evals=evaluationRuns)#100

print best