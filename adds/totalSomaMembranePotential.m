function Tsmp = totalSomaMembranePotential(t, tf, Rm, Rd, Wi, WiMat, DiMat, K, K2, stc)
% Total Soma Membrane Potential: Calculates Tsmp
%   Need some help currently on interpreting if the formulas I am using and
%   getting from Mathematica are correct.  It seems I am doing something
%   wrong and I am trying to figure out the right formulas.  In this code I
%   have tried implementing the formula I got from mathematica as shown in
%   the following link: http://i.imgur.com/4z0efXD.jpg
%   Some things in particluar I recognize are that the K and K2 variables
%   each being 'summation index of symbolic sum' seem that I am not
%   implementing them right.  I am somehow not implementing Matlab's integral 
%   function right as well.  I will try to reach out to sources to get help
%   on this but any other help I can get is greatly appreciated.
%
%   This function uses the total post synaptic current equation
%   The formula in the above link is used in Mathematica which combines formula 
%   found for solving for i (current) in the Spsc and Dpsc.
%
%   Calculates Tsmp based on t, tf, Rm, Rd, Wi, WiMat, DiMat, K, K2, stc
%   t = spike time
%   tf = pre-synaptic spike time
%   Rm = Soma resistance
%   Rd = Dendritic resistance
%   Wi = Weight attached to the dendrite
%   WiMat = Matrix of Weights
%   DiMat = Matrix of pre-synaptic spike times filtered as dirac pulses
%   size(WiMat) is used to determine the total number of input neurons
%   K = Summation index Mathematica described
%   K2 = Summation index Mathematica described
%   stc = Soma time constant

sigmaTotal = 0;
indicesInWeightMat = size(WiMat, 2);

for i=1:indicesInWeightMat
    sigmaTotal = sigmaTotal + (WiMat(i)*DiMat(i));
endfor

sigmaWiDi = (sigmaTotal / indicesInWeightMat);

K = 1;
K2 = 6290;

% Todo: acertain if K' is included in the integral.  Add in K' instead
% of K.
function outerIntegralFunct
    ((-3*Rm*t-log(tf/stc)*Rd*Rm*Wi*heaviside(-tf+K2)-Rm*t*integral(@innerIntegralFunct, 1, K2)*K)/(stc^2));
endfunction

function innerIntegralFunct
    ((log(K/stc).*sigmaWiDi)/stc);
endfunction

Tsmp = log(-t/stc)+(log(-t/stc)*integral(@outerIntegralFunct, 1, t)*K2);

endfunction
