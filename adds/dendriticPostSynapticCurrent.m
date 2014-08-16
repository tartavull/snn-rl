function Dpsc = dendriticPostSynapticCurrent(timeStep, DpscOld , Rd, Wd, Td , tf)
% Dendritic Post Synaptic Current: Calculates Dpsc
%   This function uses the total post synaptic current equation
%   Tau*Derivative[i[t], t] == - i[t] + r*w*DiracDelta[t - tf]
%   and solves for i (current).  Mathematica was used with the
%   formula DSolve[\[Tau]*D[i[t], t] == - i[t] + r*w*DiracDelta[t - tf], i, t]
%   for the result.
%
%   Calculates Dpsc based on Rd, Wi, t, tf, stc
%   Rd = Resistance dendritic
%   Wd = Weight attached to the dendrites
%   tf = pre-synaptic spike time
%   Td = Dendritic time constant

%integral from
%Dpsc = log(-t/stc+((log(-t/stc+tf/stc)*Rd*Wi*heaviside(t-tf))/stc))

Dirac = [tf tf tf tf] .* 1.8; %We are multypling it by 1.8 to get some spikes, 
% other posibility would be to increas the intial weights

Dpsc =  DpscOld + timeStep * ((-DpscOld + Rd .* Wd .* Dirac )./Td);

