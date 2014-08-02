function Dpsc = dendriticPostSynapticCurrent(Rd, Wi, t, tf, stc)
% Dendritic Post Synaptic Current: Calculates Dpsc
%   This function uses the total post synaptic current equation
%   Tau*Derivative[i[t], t] == - i[t] + r*w*DiracDelta[t - tf]
%   and solves for i (current).  Mathematica was used with the
%   formula DSolve[\[Tau]*D[i[t], t] == - i[t] + r*w*DiracDelta[t - tf], i, t]
%   for the result.
%
%   Calculates Dpsc based on Rd, Wi, t, tf, stc
%   Rd = Resistance
%   Wi = Weight attached to the dendrite
%   t = spike time.
%   tf = pre-synaptic spike time
%   stc = Soma time constant

Dpsc = log(-t/stc+((log(-t/stc+tf/stc)*Rd*Wi*heaviside(t-tf))/stc))
