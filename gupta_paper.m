close all;
clear all;
addpath('adds','datasets');
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

%Construct input Monitor
monitor = Monitor;
handle = figure(1);
monitor.setSubPlot(handle,3,2,2);
%Construct Char Monitor
charMonitor = Monitor;
charMonitor.setSubPlot(handle,3,2,1);
charMonitor.setPlotType('char');
%Construct LIK Spike Monitor
spikeMonitor = Monitor;
spikeMonitor.setSubPlot(handle,3,2,[3 4]);
spikeMonitor.setPlotType('lines3d');

%Construct adds Spike Monitor
addsMonitor = Monitor;
addsMonitor.setSubPlot(handle,3,2,[5 6]);
addsMonitor.setPlotType('lines3d');

%Initialize the net
architecture;

time = 0;
for epochIndex = 1:epochs
    for dictionaryIndex = 1:length(Dictionary)
        %Each character is presented one at a time sequentially during
        %the training process
        charCounter= mod(dictionaryIndex,length(Dictionary))+1;
        charMatrix = Dictionary{charCounter,2};
        input = reshape(charMatrix,[],1);
        
        %Record and plot
        charMonitor.record(time,charMatrix);
        charMonitor.plot();
        
        monitor.record(time,input);
        monitor.plot();
        
        for stepIndex = 0:timeStep:presentationTime          
            time = time + timeStep;
            
            %Integrate and fire layer
            likI = input * stimulusIntensity;
            likFired=find(likV >= firingThreshold);    % indices of spikes
            likV = likV + timeStep  .* 1/capacitance * (likI - likV./ resistance);
            likV(likFired) = restPotential; %Set to resting potential
            
            spikeMonitor.record(time,likV);
            spikeMonitor.plot();
            
			%Call Second layer
            LikDirac = zeros(15,1);
            LikDirac(likFired) = 1;
            
            tauDendritic = timeConstant(tauMax, tauMin , weightsDendritic);
            resistenceDendritic = resistanceComputation(tauDendritic, firingThreshold , resistenceMembrane, tauMembrane);
            
            addsFired=find(voltagesMembrane > firingThreshold);
            addsDirac = zeros(1,4);
            addsDirac(addsFired) = 1;
            %firings=[firings; addsDirac];
            
            
            currentDendritic = currentDendritic + timeStep * ((-currentDendritic + resistenceDendritic .* weightsDendritic .* [LikDirac LikDirac LikDirac LikDirac])./tauDendritic);
            currentSomatic = currentSomatic + timeStep * ((-currentSomatic + sum(weightsSomatic .* [addsDirac; addsDirac; addsDirac; addsDirac],1))./tauSomatic);
            
            voltagesMembrane =  voltagesMembrane + timeStep * ((-voltagesMembrane + resistenceMembrane .* ( sum(currentDendritic,1) + currentSomatic))/tauMembrane);
            voltagesMembrane(addsFired) = restPotential;
            
            %updateWeights
            addsMonitor.record(time,voltagesMembrane);
            addsMonitor.plot();
             
        end
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
end

