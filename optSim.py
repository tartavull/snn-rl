'''
This software uses hyperopt to optimize the parameters of the ssn-rl simulation.
Parameters are passed to the sim as arguments for running the program.
'''

from hyperopt import fmin, tpe, hp
import sys
from subprocess import * 
from decimal import Decimal

simFile = "gupta_paper_further_formulas.py"
optResultsFile = "optimizationResults.txt"
testingRunTime = 7900
minRandVals = ['minRand', 0.01, 0.50]
maxRandVals = ['maxRand', 0.51, 1.00]
posReinfVals = ['posReinfVal', 0.0, 10.0]
negReinfVals = ['negReinfVal', 0.0, 10.0]
trainTimeVals = ['trainingRunTime', 33000, 67000]
#trainTimeVals = ['trainingRunTime', 11, 33]
space = [hp.uniform(minRandVals[0],minRandVals[1],minRandVals[2]), hp.uniform(maxRandVals[0],maxRandVals[1],maxRandVals[2]), \
         hp.uniform(posReinfVals[0],posReinfVals[1],posReinfVals[2]), hp.uniform(negReinfVals[0],negReinfVals[1],negReinfVals[2]), \
         hp.uniform(trainTimeVals[0],trainTimeVals[1],trainTimeVals[2])]
evaluationRuns = 500;

# Clear and open output file
optResults = open(optResultsFile, 'w')
optResults.write("Parameters used in run:\n")
optResults.write("testingRunTime:\t\t"+str(testingRunTime)+"\tNote: testingRunTime parameter passing not currently working\n")
optResults.write(str(minRandVals)+"\n"+str(maxRandVals)+"\n"+str(posReinfVals)+"\n"+str(negReinfVals)+"\n"+str(trainTimeVals))
optResults.close()

def snnrlWithParams(kwargs):
    kwargs = kwargs.split(", ")
    print "out:\n\n",kwargs
    simRun = Popen(["python", simFile]+kwargs, stdout=PIPE)
    simRun.wait()
    results = simRun.communicate()[0].strip()
    print "\n\nresults:\n"+results+"\n\n"
    if str(results) != '': results = Decimal(results, '.1f')
    return results

def objective(args):
    minRand, maxRand, posReinfVal, negReinfVal, trainingRunTime = args
    # train
    randVal = str(minRand)+'-'+str(maxRand)
    snnrlWithParams("evaluateClassifier=False, runTime="+str(int(trainingRunTime))+", standardPrint=False, verbosePrint=False, randomization="+str(randVal)+", posReinf="+str(posReinfVal)+", negReinf="+str(negReinfVal))
    # test
    precisionPerc = snnrlWithParams("evaluateClassifier=True, optResultsFile=\'"+str(optResultsFile)+"\', runTime=\'"+str(int(trainingRunTime))+"\', standardPrint=False, verbosePrint=False, randomization="+str(randVal)+", posReinf="+str(posReinfVal)+", negReinf="+str(negReinfVal))
    precisionForOptimizing = Decimal(1, '.1f') - precisionPerc
    print "\nprecisionForOptimizing:\t",precisionForOptimizing
    return {'loss': precisionForOptimizing}

best = fmin(objective,
    space,
    algo=tpe.suggest,
    max_evals=evaluationRuns)#100

print best