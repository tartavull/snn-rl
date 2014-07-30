function Ntsmp = neuronsTSMP(charMatrix, W, O, Tmin, Tmax, stc, Rm)
% Individual Neurons Total Soma Membrane Potential: Calculates Ntsmp for Each Neuron
%   This is a function that calculates the total soma membrane potential for 
%   each individual neuron.  This code is a work in progress rough version that 
%   is not fully implemented, mainly I am showing what I am currently working 
%   on by posting this.  I understand that the Dpsc, Spsc, and Tsmp functions 
%   need to have rearranged calculations to solve for current and I am working 
%   on that.  Also, elements such as spike times (t), stmp weight adjustments, 
%   and other elements need to be added.
%   charMatrix = charactor matrix
%   W = Neuron weight matrix
%   O = Omega threshold value for spike
%   Tmin = time constant minimum
%   Tmax = time constant maximum
%   stc = Soma time constant
%   Rm = soma resistance

% todo add stmp weight adjustment code

spikeTimes = 1;
% todo: make real implementation of spike times
indicesInCharMatrix = size(charMatrix,1) * size(charMatrix,2);
DiMat = zeros(1,indicesInCharMatrix);
for i = indicesInCharMatrix
    DiMat(i) = dirac2(spikeTimes, spikeTimes);
end

%For each neuron calculate the total soma membrane potential
for Y = 1:size(charMatrix,1)
    for X = 1:size(charMatrix,2)
        charMatrixIndex = (X+(Y*size(charMatrix,2)));
        Rd = resistance(timeConstant(Tmax, Tmin, W(charCounter,charMatrixIndex)), 10, 80, 30);
        
        Dpsc = dendriticPostSynapticCurrent(charMatrix(Y,X), Rd, W(charCounter,charMatrixIndex), dirac2(spikeTimes, spikeTimes), Td);
        
        Spsc = somaticPostSynapticCurrent(W(charCounter,:), DiMat, Ts, t);
        
        Ntsmp = totalSomaMembranePotential(t, Rm, Spsc, Dpsc, Tm);
    end
end
