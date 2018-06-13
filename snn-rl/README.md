Reinforcement Learning for Spiking Neural Networks
======
Check out some of the latest results [at this link](http://nbviewer.ipython.org/github/tartavull/snn-rl/blob/master/FFSSN.ipynb) .

See the bottom of this document for file descriptions and organization of the project.

from wikipedia: "Spiking neural networks (SNNs) fall into the third generation of neural network models,
increasing the level of realism in a neural simulation. 
In addition to neuronal and synaptic state, SNNs also incorporate the concept of time into their operating model.

In practice, there is a major difference between the theoretical power of spiking neural networks and what has been demonstrated. 
They have proved useful in neuroscience, but not (yet) in engineering. As a result there has been little application of large scale 
spiking neural networks to solve computational tasks of the order and complexity that are commonly addressed 
using rate coded (second generation) neural networks. In addition it can be difficult to adapt second generation 
neural network models into real time, spiking neural networks (especially if these network algorithms are defined in 
discrete time). It is relatively easy to construct a spiking neural network model and observe its dynamics.
It is much harder to develop a model with stable behavior that computes a specific function.
"

We are trying to make the firsts steps to change that, by using a SNN to learn a very simple dataset which is compose 
by only 4 characters, described by 5 x 3 pixels.

![dataset](https://raw.githubusercontent.com/tartavull/snn-rl/master/img/readme_1.png)

We input the image as a constant current to the first layer, which contains 15 nuerons. The neurons of the first layer recive a constant current, making them fire about 3 times, while each character is shown.
The second layer contains 4 neurons, if the trainning is succesfull each neuron will be able to learn to represent one character,
spiking 3 times when the character is shown, as the image below.

![scatterplot](https://raw.githubusercontent.com/tartavull/snn-rl/master/img/readme_2.png)

======
Authors:  
* Nate Sutton
* Ignacio Tartavull

This is **work in progress** and is currently to be treated as a Proof of Concept. If you find this project interesting please **join us** and help out.

======
Getting Started

1. introduction.ipynb has some further explanations of what we are trying to achive. We are basing our work in [this paper](http://www.personal.psu.edu/lnl/papers/Gupta_Long_2007.pdf).

2. simulation.ipynb and gupta_paper_further_formulas.py are where the actual simulation is performed. We are using a clock-driven simulator for spiking neural networks called [brian](https://github.com/brian-team/brian).

3. analysis.ipynb is what we run after the simulation is finished, it helps us to analyse the performance of the network.

======
File descriptions and organization

/gupta_paper_further_formulas.py main file for project

/FFSSN.ipynb ipython notebook overviewing work

/experimentalComponents code exploring new functionality for possible use in the project

/furtherFormulas modules used in the main file

/img images

/notebooks notebooks to describe the project

/originalVersion initial version of the work, the new version is in gupta_paper_further_formulas.py 

/utilities programs used for processing files to prepare them for use in the project