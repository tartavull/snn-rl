
close all;
clear all;
addpath('adds','datasets');
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

debugScript = false;
if (debugScript)
    %Construct input Monitor
    meshMonitor = Monitor;
    handle = figure(1);
    meshMonitor.setSubPlot(handle,3,2,2);
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
end

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
        
        if (debugScript)
            %Record and plot
            charMonitor.record(time,charMatrix);
            charMonitor.plot();
            meshMonitor.record(time,input);
            meshMonitor.plot();
        end
        
        for stepIndex = 1:presentationTime*timeStep          
            time = time + timeStep;
            
            %Integrate and fire layer
            likI = input * stimulusIntensity;
            likFired=find(likV >= firingThreshold);    % indices of spikes
            likV = likV + timeStep  .* 1/capacitance * (likI - likV./ resistance);
            likV(likFired) = restPotential; %Set to resting potential
            
            if (debugScript)
                spikeMonitor.record(time,likV);
                spikeMonitor.plot();
            end
            
			%Call Second layer
            likDirac = zeros(15,1);
            likDirac(likFired) = 1;
            likFirings(:,stepIndex) = likDirac; %Something isn't working here, almost all zeros check Monitor.plotFirings
            likLastTimeFired = [likLastTimeFired likLastTimeFired(:,end)];
            likLastTimeFired(likFired,end) = time;
            
            addsFired=find(voltagesMembrane > firingThreshold);
            addsDirac = zeros(1,4);
            addsDirac(addsFired) = 1;
            addsFirings(:,stepIndex) = addsDirac;
            addsLastTimeFired = [addsLastTimeFired addsLastTimeFired(:,end)];
            addsLastTimeFired(addsFired,end) = time;
            
            tauDendritic = timeConstant(tauMax, tauMin , weightsDendritic);
            resistenceDendritic = resistanceComputation(tauDendritic, firingThreshold , resistenceMembrane, tauMembrane);        
            
            currentDendritic = currentDendritic + timeStep * ((-currentDendritic + resistenceDendritic .* weightsDendritic .* [likDirac likDirac likDirac likDirac])./tauDendritic);
            currentSomatic = currentSomatic + timeStep * ((-currentSomatic + sum(weightsSomatic .* [addsDirac; addsDirac; addsDirac; addsDirac],1))./tauSomatic);
            
            voltagesMembrane =  voltagesMembrane + timeStep * ((-voltagesMembrane + resistenceMembrane .* ( sum(currentDendritic,1) + currentSomatic))/tauMembrane);
            voltagesMembrane(addsFired) = restPotential;
                       
            if (debugScript)
                addsMonitor.record(time,voltagesMembrane);
                addsMonitor.plot();
            end
            
            %updateWeights
            %use  addsLastTimeFired -  likLastTimeFired
         
            

            
            
             
        end
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
end

