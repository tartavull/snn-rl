#references: http://stackoverflow.com/questions/3430372/how-to-get-full-path-of-current-files-directory-in-python
# http://www.tutorialspoint.com/python/python_reg_expressions.htm

import IPython.nbformat.current as nbf
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter
import os
import re

inputDirectoryPath = '/furtherFormulas/'
currentDirectory = os.getcwd()
pythonSuffix = '.py'
outputText = ''
inputDirectory =  os.listdir(currentDirectory+inputDirectoryPath)

for fileName in inputDirectory:
	if fileName.endswith(pythonSuffix) and fileName != '__init__.py':
		outputText = '';
		outputFile = (currentDirectory+inputDirectoryPath+fileName[:(len(fileName)-3)]+'.ipynb')
		f = open(currentDirectory+inputDirectoryPath+fileName,"r")
		lines = f.readlines()
		f.close()

		outputText += '# <codecell>\r\n'

		for line in lines:
			match = re.match(r'.*def (.*)\(.*:.*', line)
			if match != None and fileName != '3dAnimScatterPlotHdf5.py' and fileName != '3dBarWTauRAnim.py' and fileName != '3dBarChartAnim.py':
				outputText += '\r\n# <markdowncell>\r\n<a id=\''+match.group(1)+'\'></a>'+\
				'<div style=\'font-size:1.7em;text-decoration:underline;font-weight:bold\'>'+match.group(1)+\
				'</div>\r\n# <codecell>\r\n'
			outputText += line

		nb = nbf.reads(outputText, 'py')
		print 'out',outputFile
		nbf.write(nb, open(outputFile, 'w'), 'ipynb')

print 'done'