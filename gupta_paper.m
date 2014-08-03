close all;
clear all;
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

%Construct input Monitor
monitor = Monitor;
handle = figure(1);
monitor.setSubPlot(handle,2,2,2);
%Construct Char Monitor
charMonitor = Monitor;
charMonitor.setSubPlot(handle,2,2,1);
charMonitor.setPlotType('char');
%Construct LIK Spike Monitor
spikeMonitor = Monitor;
spikeMonitor.setSubPlot(handle,2,2,[3 4]);
spikeMonitor.setPlotType('lines3d');

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
            likFired=find(likV >= threshold);    % indices of spikes
            likV = likV + timeStep  .* 1/capacitance * (likI - likV./ resistance);
            likV(likFired) = restPotential; %Set to resting potential
            
            spikeMonitor.record(time,likV);
            spikeMonitor.plot();
            
			%Call Second layer
            addsFired=find(addsV >= threshold);
            currentDendritic = currentDendritic + timeStep * ((-currentDendritic + Rd .* weigthsDendritic .* [likFired likFired likFired])./tauDendritic);
            currentSomatic = currentSomatic + timeStep * ((-currentSomatic + sum(weightsSomatic .* [addsFired addsFired addsFired addsFired],2))./tauSomatic);
            
            
        end
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
end

