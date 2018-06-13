#references: http://stackoverflow.com/questions/3430372/how-to-get-full-path-of-current-files-directory-in-python
# http://www.tutorialspoint.com/python/python_reg_expressions.htm

import IPython.nbformat.current as nbf
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter
import os
import re

currentDirectory = os.getcwd()
fileName='OfflineFFSSN.ipynb'
outputFile='FFSSN.ipynb'

outputText = '';
f = open(currentDirectory+'/'+fileName,"r")
lines = f.readlines()
f.close()

for line in lines:
	match = re.match(r'(.*)(notebooks)(.*)', line)
	if match != None:
		outputText += match.group(1)+'github/tartavull/snn-rl/blob/master'+match.group(3)
	else:
		outputText += line

#nb = nbf.reads(outputText, 'ipynb')
#nbf.write(nb, open(outputFile, 'w'), 'ipynb')
f = open(currentDirectory+'/'+outputFile,"w")
f.write(outputText)
f.close()

print 'done'