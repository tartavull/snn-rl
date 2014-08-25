function Wnew = newWeight(Dw, Wold, Wmin, Wmax, n)
% New Weight: Trains weights based on parameters such as
%   the learning rate.
%   Dw = delta weight function output.
%   Wold = origional weight
%   Wmin = minimum weight parameter
%   Wmax = maximum weight parameter
%   n = learning rate

positives = find(Dw >= 0);
Wnew(positives) = Wold(positives) + n .* Dw(positives) .* ( Wmax - Wold(positives) );

negatives = find(Dw < 0);
Wnew(negatives) = Wold(negatives) + n .* Dw(negatives) .* ( Wold(negatives) - Wmin );
endfunction
