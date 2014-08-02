function Wnew = newWeight(Dw, Wold, Wmin, Wmax, n)
% New Weight: Trains weights based on parameters such as
%   the learning rate.
%   Dw = delta weight function output.
%   Wold = origional weight
%   Wmin = minimum weight parameter
%   Wmax = maximum weight parameter
%   n = learning rate

if (Dw >= 0)
   Wnew = Wold+n*Dw*(Wmax-Wold);
elseif (Dw < 0)
   Wnew = Wold+n*Dw*(Wold-Wmin);
end
