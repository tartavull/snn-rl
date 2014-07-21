% daspnet.m: Spiking network with DA-modulated STDP. Based on spnet.m
% Created by Eugene M.Izhikevich.                November 3, 2004
%
% This program reproduces the experiment in Fig.1 in
% Izhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage 
% of STDP and Dopamine Signaling. Cerebral Cortex, 10.1093/cercor/bhl152 
%
% n1 - the presynaptic neuron. syn is the synapse to be reinforced.
% Plot: top - spike raster. Bottom left - synaptic strength (blue), the
% eligibility trace (green), and the rewards (red x). Bottom right - the
% distribution of synaptic weights with the chosen synapse marked by red dot.

M=100;                 % number of synapses per neuron
D=1;                   % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];
d=[   8*ones(Ne,1);    2*ones(Ni,1)];
sm=4;                 % maximal synaptic strength

post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
s=[ones(Ne,M);-ones(Ni,M)];         % synaptic weights
sd=zeros(N,M);                      % their derivatives
for i=1:N
  if i<=Ne
    for j=1:D
      delays{i,j}=M/D*(j-1)+(1:M/D);
    end;
  else
    delays{i,1}=1:M;
  end;
  pre{i}=find(post==i&s>0);             % pre excitatory neurons
  aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
end;
STDP = zeros(N,1001+D);
v = -65*ones(N,1);                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

%---------------
% new stuff related to DA-STDP
T=3600;         % the duration of experiment
DA=0;           % level of dopamine above the baseline
rew=[];

n1=1;           % presynaptic neuron
syn=1;          % the synapse number to the postsynaptic neuron
n2=post(n1,syn) % postsynaptic neuron
s(n1,syn)=0;    % start with 0 value

interval = 20;  % the coincidence interval for n1 and n2
n1f=-100;       % the last spike of n1
n2f=[];         % the last spike of n2
shist=zeros(1000*T,2);
%--------------

for sec=1:T                             % simulation of 1 day
  for t=1:1000                          % simulation of 1 sec
    I=13*(rand(N,1)-0.5);               % random thalamic input 
    fired = find(v>=30);                % indices of fired neurons
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=0.1;
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-1.5*STDP(ind,t+D)';
      k=k-1;
    end;
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % for numerical 
    v=v+0.5*((0.04*v+5).*v+140-u+I);    % stability time 
    u=u+a.*(0.2*v-u);                   % step is 0.5 ms
    STDP(:,t+D+1)=0.95*STDP(:,t+D);     % tau = 20 ms
   
    DA=DA*0.995;
    if (mod(t,10)==0)
        s(1:Ne,:)=max(0,min(sm,s(1:Ne,:)+(0.002+DA)*sd(1:Ne,:)));
        sd=0.99*sd;
    end;
    if any(fired==n1)
        n1f=[n1f,sec*1000+t];
    end
    if any(fired==n2)
        n2f=[n2f,sec*1000+t];
        if (sec*1000+t-n1f(end)<interval) & (n2f(end)>n1f(end))
            rew=[rew,sec*1000+t+1000+ceil(2000*rand)];
        end;
    end     
    if any(rew==sec*1000+t)
        DA=DA+0.5;
    end;
    shist(sec*1000+t,:)=[s(n1,syn),sd(n1,syn)];

  end;
% ---- plot -------
  subplot(2,1,1)
  plot(firings(:,1),firings(:,2),'.');
  axis([0 1000 0 N]); 
  subplot(2,2,3);
  plot(0.001*(1:(sec*1000+t)),shist(1:sec*1000+t,:), 0.001*rew,0*rew,'rx');
  subplot(2,2,4);
  hist(s(find(s>0)),sm*(0.01:0.01:1)); % only excitatory synapses
  hold on; plot(s(n1,syn),0,'r.'); hold off;
  drawnow;
% ---- end plot ------
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
end;
