%Simulation time parameters
epochs = 10; %an epoch means one full presentation % of all the characters 100 
presentationTime = 300; %Each character is presented for 300ms
timeStep = 0.2; %time step for the simulation is 0.2 ms

%LIF parameters
restPotential = 0;
likV = restPotential * ones(15,1);  % Initial values of voltage
stimulusIntensity = 20;
capacitance = 20;
resistance = 0.8;

%Weights
weightsDendritic = (0.5 * rand(1) + 0.5) * ones(15,4);
weightsSomatic = -(0.5 * rand(1) + 0.5) * ones(4,4); % We dont know if this initialization is correct
weightMaxExcitatory = 1;
weightMinExcitatory = 0;
weightMaxInhibitory = 0;
weightMinInhibitory = -1;
aPLus = 0.1;
aMinus = -0.105;

%Time constants
tauMax = 30;
tauMin = 2;
tauDendritic = zeros(15,4);
tauSomatic = 2;
tauMembrane = 2;
tauPlus = 1;
tauMinus = 1;

%Currents
currentDendritic = zeros(15,4);
currentSomatic = zeros(1,4);

%Firings
LikFirings  = zeros(length(Dictionary) * epochs * presentationTime / timeStep , 15);
addsFirings = zeros(length(Dictionary) * epochs * presentationTime / timeStep , 4);

%Resistences
resistenceDendritic = zeros(15,4);
resistenceMembrane = 80;

%Voltages
voltagesMembrane = zeros(1,4);

%Firing threshold
firingThreshold = 10;

%Learning rate
learningRate = 0.1;
learningRateDecay = 0.05;