% LIF MATLAB model - based on code
% and lecture slides of Prof. Kevin Gurney

%% CLEANUP
clear all
clc
clf

%% PARAMETER SETUP
% membrane constants
tau = 0.020;
R = 3e7;
% resting potential
E = -0.070;
% threshold for a spike
theta = -0.030;
% change in time - window
dt = 0.0001; %should be less than tau
% total milliseconds to run for
T = 0.16;
% total number of steps
no_steps = round(T ./ dt);
% time - for plotting
time = linspace(0, T, no_steps + 1);
% external events array
psc = [];
 % external inhibition event
ipsc = [];
% set up PSC events - external current give rise to inc v.  
psc(1) = 0.05;
psc(2) = 0.09;
ipsc(1) = 0.04;
% number of external events
no_pcs = length(psc);
no_ipcs = length(ipsc);
tau_s = 0.003; %should be less than tau
Q = 40e-12;
I_0 = Q ./ tau_s;
% index for each ext. event? IE psc(1)/dt, psc(2)/dt...
index_pscs = round(psc ./ dt);
index_ipscs = round(ipsc ./ dt);

% voltage matrix
V = zeros(1, no_steps + 1);
% initial voltage of membrane - put at rest
V(1) = -0.07;
% current
I = zeros(1, no_steps + 1);
% initial current
I(1) = 0;
% noise - sampling from norm/Gaussian dist, 1*no_steps matrix
% random(A,B,C,[D,E])- what does B do!?
randI = 3e-9 .* random('Normal', 0, 1.5, [1, no_steps]);
% time since last spike
t_spike = 0;
% absolute refactory period
arp = 0.01; %this being 0.02 doesn't work.
% total number of spikes
no_spikes = 0;

%% SIMULATION
for i=1:no_steps
    for k=1:no_pcs
        % if step sits on an external event step...
        if i == index_pscs(k)
            % increase the current
            I(i) = I(i) + I_0;
        end
    end
    for x=1:no_ipcs
        if i == index_ipscs(x)
            I(i) = I(i) - I_0;
        end
    end
    I(i+1) = I(i)-(dt/tau_s).*I(i);
    
    dV =(dt/tau).*(E-V(i)+I(i).*R + randI(i).*R);
    % update without noise
    % dV =(dt/tau).*(E-V(i)+I(i).*R);
    V(i+1) = V(i) + dV;
    % spike
    if (V(i+1) > theta) 
        if no_spikes>0
            % check we're not in absolute refac period
            if (time(i)>=(t_spike+arp))
                % reset voltage
                V(i+1) = E;
                % Eleni's trick of making spikes look nice
                V(i) = -0;
                % record spike
                t_spike = time(i);
                % increment spike count
                no_spikes = no_spikes+1;
            end
        else
            % no spikes yet - no need to check
            V(i+1) = E;
            % Eleni's trick of making spikes look nice
            V(i) = -0;
            % record spike
            t_spike = time(i);
            % increment spike count
            no_spikes = no_spikes+1;
        end
    end
end
%% PLOTTING
figure(1)
plot(time, V, 'color', 'black');
hold on
plot(time, theta, '--r');
line1 = 'Leaky integrate-and-fire model with noise';
line2 = sprintf('theta: %.3fV, refactory period: %.2fs',theta, arp);
title({line1, line2});
xlabel('time (s)')
ylabel('voltage (V)')
grid on