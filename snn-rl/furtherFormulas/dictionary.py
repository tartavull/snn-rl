'''
Updated 03/24/18: matplotlib imports added to avoid error with newer matplotlib versions. Also fixed 
dictionary for newer numpy version issues.
'''

import matplotlib
matplotlib.use('qt4agg') 
from matplotlib.pyplot import *#as plt
from brian2 import *

#In the matlab implementation the fisrt layer consited of 15 leaky integrate and fire in which 
#the input was a constant current in case the pixel was "on" or a null current in case the pixel 
#was "off". In the case the pixel was on, the nueron that represented that pixel would fire around 3 times 
#in that epoch.
#In this implementation, we use the class SpikeGeneratorGroup, to generate spikes at spicified times. 
#By doing it this way, we have more control over the real neurons ('second layer') we are trying to train.
#and a simpler implementation. This "dictionary" class is responsable for generating the array which contains
#the firing time for the 15 neurons.
 
class dictionary():
	def spikeTimes(self, dictionaryLongitude, spikeInterval, spikesPerChar, epochs):

		dictionary = [['A',np.array([[1, 0, 1],[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0]])],
					    ['B',np.array([[0, 0, 1],[0, 1, 0],[0, 0, 1],[0, 1, 0],[0, 0, 1]])],
					    ['C',np.array([[0, 0, 0],[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 0, 0]])],
					    ['D',np.array([[0, 0, 1],[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 1]])],
					    ['E',np.array([[0, 0, 0],[0, 1, 1],[0, 0, 1],[0, 1, 1],[0, 0, 0]])],
					    ['F',np.array([[0, 0, 0],[0, 1, 1],[0, 0, 1],[0, 1, 1],[0, 1, 1]])],
					    ['G',np.array([[0, 0, 0],[0, 1, 1],[0, 1, 1],[0, 1, 0],[0, 0, 0]])],
					    ['H',np.array([[0, 1, 0],[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0]])],
					    ['I',np.array([[0, 0, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1],[0, 0, 0]])],
					    ['J',np.array([[0, 0, 0],[1, 1, 0],[1, 1, 0],[0, 1, 0],[0, 0, 0]])],
					    ['K',np.array([[0, 1, 0],[0, 0, 1],[0, 1, 1],[0, 0, 1],[0, 1, 0]])],
					    ['L',np.array([[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 0, 0]])],
					    ['M',np.array([[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0]])],
					    ['N',np.array([[0, 1, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 1, 0]])],
					    ['O',np.array([[1, 0, 1],[0, 1, 0],[0, 1, 0],[0, 1, 0],[1, 0, 1]])],
					    ['P',np.array([[0, 0, 0],[0, 1, 0],[0, 0, 0],[0, 1, 1],[0, 1, 1]])],
					    ['Q',np.array([[1, 0, 1],[0, 1, 0],[0, 1, 0],[1, 0, 1],[1, 1, 0]])],
					    ['R',np.array([[0, 0, 1],[0, 1, 0],[0, 0, 0],[0, 0, 1],[0, 1, 0]])],
					    ['S',np.array([[1, 0, 0],[0, 1, 1],[1, 0, 1],[1, 1, 0],[0, 0, 1]])],
					    ['T',np.array([[0, 0, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1],[1, 0, 1]])],
					    ['U',np.array([[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 0]])],
					    ['V',np.array([[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0],[1, 0, 1]])],
					    ['W',np.array([[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 0],[0, 1, 0]])],
					    ['X',np.array([[0, 1, 0],[0, 1, 0],[1, 0, 1],[0, 1, 0],[0, 1, 0]])],
					    ['Y',np.array([[0, 1, 0],[0, 1, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1]])],
					    ['Z',np.array([[0, 0, 0],[1, 1, 0],[1, 0, 1],[0, 1, 1],[0, 0, 0]])]]
		dictionary = np.array(dictionary, dtype=object)
		
		self.dictionary = dictionary

		spikeArray = array([]).reshape(0,2)

		time = 0 * ms

		for indexEpoch in range (0, epochs):
			for indexDictionary in range(0,dictionaryLongitude):
				
				vector = dictionary[indexDictionary][1].reshape([1,15])
				for charSpikeIndex in range(0,spikesPerChar):

					time = time + spikeInterval
					for indexPixel in range(0,15):  
						
						if vector[0][indexPixel] == 1:
						
							spike = array([indexPixel, time])
							spikeArray = vstack([spikeArray,spike])
		
		return spikeArray 


				    
