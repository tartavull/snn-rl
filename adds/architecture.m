
%Simulation time parameters
epochs = 100; %an epoch means one full presentation % of all the characters 100 
presentationTime = 300; %Each character is presented for 300ms
timeStep = 0.2; %time step for the simulation is 0.2 ms

makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

lenAdds = length(Dictionary);
lenLik = length(reshape(Dictionary{1,2},[],1));

%LIF parameters
restPotential = 0;
likV = restPotential * ones(lenLik,1);  % Initial values of voltage
stimulusIntensity = 20;
capacitance = 20;
resistance = 0.8;

%Weights
weightsDendritic = 0.5 .* rand(lenLik,lenAdds) + 0.5;
%weightsDendritic = rand(lenLik,lenAdds);

weightsSomatic = 0.5 * rand(lenAdds,lenAdds) -0.75; 
%We should keep the diagonal as zeroes so we dont autoinhibate the spikes
weightsSomatic(1:lenAdds+1:lenAdds*lenAdds) = 0;

weightMaxExcitatory = 1;
weightMinExcitatory = 0;
weightMaxInhibitory = 0;
weightMinInhibitory = -1;
aPLus = 0.1;
aMinus = -0.105;

%Time constants
tauMax = 30;
tauMin = 2;
tauDendritic = zeros(lenLik,lenAdds);
tauSomatic = 2;
tauMembrane = 2;
tauPlus = 1;
tauMinus = 1;

%Currents
currentDendritic = zeros(lenLik,lenAdds);
currentSomatic = zeros(1,lenAdds);

%Firings
likFirings  = zeros(lenLik,length(Dictionary) * epochs * presentationTime / timeStep );
addsFirings = zeros(4,length(Dictionary) * epochs * presentationTime / timeStep);

likLastTimeFired= zeros(lenLik,1) * NaN;
addsLastTimeFired= zeros(lenAdds,1) * NaN;


%Resistences
resistenceDendritic = zeros(lenLik,lenAdds);
resistenceMembrane = 80;

%Voltages
voltagesMembrane = zeros(1,lenAdds);

%Firing threshold
firingThreshold = 10;

%Learning rate
learningRate = 0.1;
learningRateDecay = 0.05;


