function Dw = deltaWeight(Dt)
% Delta Weight: Calculates changes in weights based on the Spike Time
%   Dependent Plasticity (STDP) learning rule.
%   Calculates delta weight based on Dt
%   Dt = delta spike time.  
%   A ln constant and not function was used because it was found the ln
%   function is a non-standard matlab function and not everyone may have
%   it.

APlus = 0.1; % parameter based on appendix I
AMinus = -0.105; % parameter based on appendix I
TimeConstantPlus = 1.0; % time constant in ms parameter based on appendix I
TimeConstantMinus = 1.0; % time constant in ms parameter based on appendix I
ln = 2.71828182845904523536028747135266249775724709369995; % natural log
   
   Dw = zeros(length(Dt),1);
   negatives =  find(Dt < 0);
   positives =  find(Dt > 0);
   Dw(negatives) =  APlus .* (ln .^ (Dt(negatives)./TimeConstantPlus));
   Dw(positives) =  AMinus .* (ln .^ (-Dt(positives)./TimeConstantMinus));
end
