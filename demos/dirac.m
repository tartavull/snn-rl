function Di = dirac(Dt, t)
%Dirac function: Calculates Di
%   Uses the Dirac function as described in the book 'Theoretical 
%   Neuroscience' on page 405.  The function can be defined in other 
%   ways but this is one way to define it.  This expresses the function
%   as the limit of a square pulse.
%   Dt = delta spike time.  
%   t = spike time.
%   The formula below represents the limit of the reciprocal of Dt as
%   Dt approaches 0.  
%   Note: this code should work but it may be a couple of days before
%   I have access to the limit function to fully test it.

if (t > -Dt / 2) && (t < Dt / 2)
    Di = limit(1 / Dt, Dt, 0);
else
    Di = 0;
end
