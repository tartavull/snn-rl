function Dpsc = dendriticPostSynapticCurrent(Id, Rd, Wi, Di)
% Dendritic Post Synaptic Current: Calculates Dpsc
%   Calculates Dpsc based on Id, Rd, Wi, Di
%   Id = Dendritic current input
%   Rd = Resistance
%   Wi = Weight attached to the dendrite
%   Di = Dirac pulses calculated using the dirac function

Dpsc = -Id + Rd * Wi * Di;
