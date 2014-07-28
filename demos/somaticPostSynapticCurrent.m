function Spsc = somaticPostSynapticCurrent(Is, WiMat, DiMat, NTotal)
% Somatic Post Synaptic Current: Calculates Spsc
%   Calculates Spsc based on Is, WiMat, DiMat, NTotal
%   I did not include the matlab sigma function because some contributers
%   may not have the Symbolic Toolbox and I replaced it with a loop.
%   Is = Somatic current input
%   WiMat = Matrix of Weights
%   DiMat = Matrix of pre-synaptic spike times filtered as dirac pulses
%   NTotal is the total number of input neurons

sigmaTotal = 0;

for i=1:NTotal
    sigmaTotal = sigmaTotal + (WiMat(i)*DiMat(1));
end

sigmaWiDi = (sigmaTotal / NTotal);

Spsc = -Is + sigmaWiDi;
