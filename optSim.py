from hyperopt import fmin, tpe, hp
import gupta_paper_further_formulas as snnrl
import sys

def snnrlWithParams(*args, **kwargs):
	sys.argv = ['snnrl.py', args, kwargs]
	return snnrl.main()	

def objective(args):
    minRand, maxRand, posReinfVal, negReinfVal = args
    # train
    randVal = str(minRand)+'-'+str(maxRand)
    snnrlWithParams(evaluateClassifier=False, randomization=randVal, posReinf=posReinfVal, negReinf=negReinfVal)
    # test
    accuracyPerc = 1 - (snnrlWithParams(evaluateClassifier=True, posReinf=posReinfVal, negReinf=negReinfVal))
    print "\naccPer:\t",accuracyPerc

    return {'loss': accuracyPerc}

space = [hp.uniform('minRand', 0.6, 0.7), hp.uniform('maxRand', 0.9, 0.9), hp.uniform('posReinfVal', 1.0, 5.0), hp.uniform('negReinfVal', 1.0, 5.0)]

best = fmin(objective,
    space,
    algo=tpe.suggest,
    max_evals=3)#100

print best