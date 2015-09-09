'''
Profiling sim to find ways to quicken run time
references: 
https://zapier.com/engineering/profiling-python-boss/
http://java2s.com/Code/Python/Class/Classdestructordel.htm
'''
from hyperopt import fmin, tpe, hp
import gupta_paper_further_formulas as snnrl
import sys
import cProfile

def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func

def snnrlWithParams(*args, **kwargs):
    sys.argv = ['snnrl.py', args, kwargs]
    '''gupta_paper()
    accuracyPerc = run_gupta_paper.run_model()'''
    return snnrl.gupta_paper() 

@do_cprofile
def objective(args):
    minRand, maxRand, posReinfVal, negReinfVal = args
    # train
    randVal = str(minRand)+'-'+str(maxRand)
    trainModel = snnrlWithParams(evaluateClassifier=False, randomization=randVal, posReinf=posReinfVal, negReinf=negReinfVal)
    trainModel.run_model()
    del trainModel
    # test
    testModel = snnrlWithParams(evaluateClassifier=True, posReinf=posReinfVal, negReinf=negReinfVal)
    origAccPerc = testModel.run_model()
    accuracyPerc = 1 - (origAccPerc)
    print "\naccPer:\t",accuracyPerc
    del testModel

    return {'loss': accuracyPerc}

#objective([0.6, 0.9, 1.5, 3.0])    

space = [hp.uniform('minRand', 0.6, 0.7), hp.uniform('maxRand', 0.9, 0.9), hp.uniform('posReinfVal', 1.0, 5.0), hp.uniform('negReinfVal', 1.0, 5.0)]

best = fmin(objective,
    space,
    algo=tpe.suggest,
    max_evals=3)#100

print best