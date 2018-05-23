from brian import *

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

		dictionary = zeros(26, dtype='c1, (5,3)u1')
		dictionary[:] =[('A',[[1, 0, 1],[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0]]),
					    ('B',[[0, 0, 1],[0, 1, 0],[0, 0, 1],[0, 1, 0],[0, 0, 1]]),
					    ('C',[[0, 0, 0],[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 0, 0]]),
					    ('D',[[0, 0, 1],[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 1]]),
					    ('E',[[0, 0, 0],[0, 1, 1],[0, 0, 1],[0, 1, 1],[0, 0, 0]]),
					    ('F',[[0, 0, 0],[0, 1, 1],[0, 0, 1],[0, 1, 1],[0, 1, 1]]),
					    ('G',[[0, 0, 0],[0, 1, 1],[0, 1, 1],[0, 1, 0],[0, 0, 0]]),
					    ('H',[[0, 1, 0],[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0]]),
					    ('I',[[0, 0, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1],[0, 0, 0]]),
					    ('J',[[0, 0, 0],[1, 1, 0],[1, 1, 0],[0, 1, 0],[0, 0, 0]]),
					    ('K',[[0, 1, 0],[0, 0, 1],[0, 1, 1],[0, 0, 1],[0, 1, 0]]),
					    ('L',[[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 1, 1],[0, 0, 0]]),
					    ('M',[[0, 1, 0],[0, 0, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0]]),
					    ('N',[[0, 1, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 1, 0]]),
					    ('O',[[1, 0, 1],[0, 1, 0],[0, 1, 0],[0, 1, 0],[1, 0, 1]]),
					    ('P',[[0, 0, 0],[0, 1, 0],[0, 0, 0],[0, 1, 1],[0, 1, 1]]),
					    ('Q',[[1, 0, 1],[0, 1, 0],[0, 1, 0],[1, 0, 1],[1, 1, 0]]),
					    ('R',[[0, 0, 1],[0, 1, 0],[0, 0, 0],[0, 0, 1],[0, 1, 0]]),
					    ('S',[[1, 0, 0],[0, 1, 1],[1, 0, 1],[1, 1, 0],[0, 0, 1]]),
					    ('T',[[0, 0, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1],[1, 0, 1]]),
					    ('U',[[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 0]]),
					    ('V',[[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 1, 0],[1, 0, 1]]),
					    ('W',[[0, 1, 0],[0, 1, 0],[0, 1, 0],[0, 0, 0],[0, 1, 0]]),
					    ('X',[[0, 1, 0],[0, 1, 0],[1, 0, 1],[0, 1, 0],[0, 1, 0]]),
					    ('Y',[[0, 1, 0],[0, 1, 0],[1, 0, 1],[1, 0, 1],[1, 0, 1]]),
					    ('Z',[[0, 0, 0],[1, 1, 0],[1, 0, 1],[0, 1, 1],[0, 0, 0]])]
		
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


				    
