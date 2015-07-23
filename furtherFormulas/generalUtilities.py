#from __future__ import print_function
import numpy as np
from decimal import Decimal
import re
#import sys

def parseArgs(self, systemArgs, dictionaryLongitude):
	'''
	Aruments passed to the program are parsed and variable values are
	extracted from them.

	The join(re.findall... one liner regex command was from:
	http://stackoverflow.com/questions/23366848/getting-captured-group-in-one-line
	It parses text from within single quotes that exists in the kwargs

	'''

	#print '\n', systemArgs, '\n'
	# process command line args
	for allArgs in systemArgs:
		argsGroup = str(allArgs).split(',')
		for args in argsGroup:
			keyAndVal = args.split('=')
			# check for kwarg not arg
			if (len(keyAndVal) > 1):
				# Preprocessing
				currentKWArgName = ''.join(re.findall(r'\'(.*)\'', keyAndVal[0]))
				currentKWArgValue = ''.join(re.findall(r'\'(.*)\'', keyAndVal[1]))
				if currentKWArgName == '': currentKWArgName = keyAndVal[0]
				if currentKWArgValue == '': currentKWArgValue = keyAndVal[1]
				currentKWArgValue = currentKWArgValue.replace("}", "")
				currentKWArgValue = currentKWArgValue.replace("{", "")
				if currentKWArgValue == 'True': currentKWArgValue = True
				elif currentKWArgValue == 'False': currentKWArgValue = False						
				# Saving vals
				if currentKWArgName == "standardPrint":
					self.standardPrint = currentKWArgValue
				elif currentKWArgName == "verbosePrint":
					self.verbosePrint = currentKWArgValue	
				elif currentKWArgName == "evaluateClassifier":
					self.evaluateClassifier = currentKWArgValue
					# By default accelerate
					if currentKWArgValue == False: self.accelerateTraining = True
				elif currentKWArgName == "accelerateTraining":
					self.accelerateTraining = currentKWArgValue
				elif currentKWArgName == "randomization":
					self.minWeightsRand = Decimal(currentKWArgValue.split("-")[0], '.1f')
					self.maxWeightsRand = Decimal(currentKWArgValue.split("-")[1], '.1f')
					#print("\n+++rand found+++\t"+str(self.minWeightsRand)+"\t"+str(self.maxWeightsRand)+"\n")
					self.W = np.random.uniform(self.minWeightsRand,self.maxWeightsRand,[dictionaryLongitude,15]) # Initial weights
					W = self.W
				elif currentKWArgName == "posReinf":
					self.positiveWeightReinforcement = float(currentKWArgValue)#Decimal(currentKWArgValue, '.1f')
				elif currentKWArgName == "negReinf":
					self.negativeWeightReinforcement = float(currentKWArgValue)#Decimal(currentKWArgValue, '.1f')			
				elif currentKWArgName == "showPlot":
					self.showPlot = currentKWArgValue
				elif currentKWArgName == "runTime" or currentKWArgName == "testingRunTime":
					from brian2 import ms
					if currentKWArgName == "runTime":
						self.runTime = int(currentKWArgValue)*ms
					else:
						# TODO: testingRunTime argument passing is not working and it can be fixed later
						self.testingRunTime = int(currentKWArgValue)*ms		
				elif currentKWArgName == "optResultsFile":
					self.optResultsFile = currentKWArgValue

				if self.verbosePrint == True:
					print("\n1:  ",keyAndVal[0]," & ",currentKWArgName," 2:  ",keyAndVal[1]," & ",currentKWArgValue,"\n")

	if self.standardPrint == True:
		print("\n\nParameters that can be used when running the program are: \n\
				<trainOrTest>: if training or testing should be run\n\
				<randomization>: the amount of randomization to initialize weights for\n\
				training with.  Default value: \"randomization = 0.50-1.00\"\n\
				<posReinf>: strength of positive reinforcement.  Default value: \"1.0\"\n\
				<negReinf>: strength of negative reinforcement.  Default value: \"2.0\"\n\
				<standardPrint>: print regular output.  Note: accuracy printed irregardless\n\
				<verbosePrint>: print extra details\n\
				")

	if self.verbosePrint == True:
		print('initial Weights\n',self.W)

	#if (self.standardPrint == False and self.verbosePrint == False):
	#	print("Starting sim run for total time: "+str(self.runTime * self.runTimeScaling), file=sys.stderr)

	if self.evaluateClassifier == True: 
		if (self.optResultsFile != ""):
			optResults = open(self.optResultsFile, 'a')
			optResults.write("\nTestTime: "+str(self.runTime))
			optResults.close()

	return self	