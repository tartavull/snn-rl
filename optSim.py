from hyperopt import fmin, tpe, hp
#import gupta_paper_further_formulas as snnrl
import sys
from subprocess import * 
from decimal import Decimal

simFile = "gupta_paper_further_formulas.py"
trainingRunTime = 1100

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
    results = snnrlWithParams("evaluateClassifier=True, standardPrint=False, verbosePrint=False, posReinf="+str(posReinfVal)+", negReinf="+str(negReinfVal))
    accuracyPerc = Decimal(1, '.1f') - results
    print "\naccPer:\t",accuracyPerc

    return {'loss': accuracyPerc}

space = [hp.uniform('minRand', 0.6, 0.7), hp.uniform('maxRand', 0.9, 0.9), hp.uniform('posReinfVal', 1.0, 5.0), hp.uniform('negReinfVal', 1.0, 5.0)]

best = fmin(objective,
    space,
    algo=tpe.suggest,
    max_evals=3)#100

print best