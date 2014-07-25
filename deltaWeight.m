function Dw = deltaWeight(Dt)
% Delta Weight: Calculates changes in weights based on the Spike Time
%   Dependent Plasticity (STDP) learning rule.
%   Calculates delta weight based on Dt
%   Dt = delta spike time.  

APlus = 0.1 % parameter based on appendix I
AMinus = -0.105 % parameter based on appendix I
TimeConstantPlus = 1 % time constant in ms parameter based on appendix I
TimeConstantMinus = 1 % time constant in ms parameter based on appendix I

if (Dt < 0)
   Dw =  APlus * ln(Dt/TimeConstantPlus);
elseif (Dt > 0)
   Dw =  AMinus * ln(-Dt/TimeConstantMinus);
end
