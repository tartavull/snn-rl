'''

References:
http://stackoverflow.com/questions/89228/calling-an-external-command-in-python
https://docs.python.org/2/library/subprocess.html#module-subprocess
'''

#import subprocess as subprocess
from subprocess import * 
from decimal import Decimal

#data = call(["ls", "-l"])

#evaluateClassifier=False, randomization=randVal, posReinf=posReinfVal, negReinf=negReinfVal

#output = subprocess.check_output(["python /home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py", "evaluateClassifier=True randomization='0.5-1.0' posReinf='2.0' negReinf='3.0'"])
#output = subprocess.check_output(["python", "/home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py evaluateClassifier=False"])
#output = subprocess.check_output(["python /home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py"])
#output = subprocess.call(["python", "/home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py", "evaluateClassifier=True", "standardPrint = False", shell=True)
#output = subprocess.call("exit 1", shell=True)
#output = subprocess.call(["python", "/home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py", "evaluateClassifier=True", "standardPrint=False", "verbosePrint=False"])
testRun = Popen(["python", "/home/nmsutton/Documents/Software/SNNRL/snn-rl/gupta_paper_further_formulas.py", "evaluateClassifier=True", "standardPrint=False", "verbosePrint=False"], stdout=PIPE)
#testRun = Popen(["ls", "-l"], stdout=PIPE)
#p = Popen(["/bin/ls", "-l"], stdout=PIPE)
#(child_stdin, child_stdout) = (p.stdin, p.stdout)


#print len(data)
#print '\n\nData:\n\n\'

testRun.wait()
print '\nstout:\n'
results = testRun.communicate()[0].strip()
print results, '\n'
print Decimal(results, '.1f')*Decimal(3.0, '.1f')
#print output
