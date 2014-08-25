function Rd = resistanceComputation(Td, theta, Rm, Tm)
%resistance: Calculates resistance
%   Calculates resistance based on Time constant (Td), Firing threshold
%   (theta), Soma resistance (Rm), Soma time constant (Tm).
Tm = ones(size(Td)) .* Tm;

Rd =  (Td .* theta ./ Rm) .* (Tm ./ Td).^(Tm ./ (Tm-Td));

endfunction

