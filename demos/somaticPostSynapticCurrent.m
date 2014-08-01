function Spsc = somaticPostSynapticCurrent(WiMat, DiMat, stc, t)
% Somatic Post Synaptic Current: Calculates Spsc
%   Need some help currently on interpreting how to code 'K' from the
%   following Mathematica 'out' formula: http://i.imgur.com/Tfa14od.jpg
%   appears to be this function http://mathworld.wolfram.com/K-Function.html
%   Trying to figure out how can that be coded in Matlab?
%   I am also still working on getting the integral function to work.
%
%   This function uses the total post synaptic current equation
%   Tau*D[i[t], t] == - i[t] + sum[w*DiracDelta[t - tf]]
%   and solves for i (current).  Mathematica was used with the
%   formula DSolve[\[Tau]*D[i[t], t] == - i[t] + sum[w*DiracDelta[t - tf]], i, t]
%
%   Calculates Spsc based on WiMat, DiMat, stc, t
%   I did not include the matlab sigma function because some contributers
%   may not have the Symbolic Toolbox and I replaced it with a loop.
%   WiMat = Matrix of Weights
%   DiMat = Matrix of pre-synaptic spike times filtered as dirac pulses
%   size(WiMat) is used to determine the total number of input neurons

sigmaTotal = 0;
indicesInWeightMat = size(WiMat, 2);

for i=1:indicesInWeightMat
    sigmaTotal = sigmaTotal + (WiMat(i)*DiMat(i));
end

sigmaWiDi = (sigmaTotal / indicesInWeightMat);

K = 1;

% Todo: acertain if K' is included in the integral.  Add in K' instead
% of K.
function integralFunct
    disp(((log(1/stc).*sigmaWiDi)/stc).*K);
end

Spsc = log(-t/stc)+(log(-t/stc)*integral(@integralFunct, 1, t));

end
