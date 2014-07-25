function Di = dirac(Dt, t)
%Dirac function: Calculates Di
%   Uses the Dirac function as described in the book 'Theoretical 
%   Neuroscience' on page 405.  The function can be defined in other 
%   ways but this is one way to define it.  This expresses the function
%   as the limit of a square pulse.
%   Dt = delta spike time.  
%   t = spike time.
%   The origional dirac formula represents the limit of the reciprocal 
%   of Dt as Dt approaches 0.  The limit function is only availible in the
%   Symbolic Math Toolbox and some other may not have that so I have
%   excluded it here.
%   TODO: Now looking into another way to create the limit function or what
%   effect removing it may cause.

if (t > -Dt / 2) && (t < Dt / 2) && Dt ~= 0
    Di = 1 / Dt;
else
    Di = 0;
end
