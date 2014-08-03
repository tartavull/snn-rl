close all;
clear all;
addpath('adds','datasets');
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


%Leaky integrate and fire intialization
restPotential = 0;
likV = restPotential * ones(15,1);  % Initial values of voltage
stimulusIntensity = 20;
threshold = 10;
capacitance = 20;
resistance = 0.8;

% Second layer intialization
outputNeurons = 4;
voltageSecLayer = zeros(outputNeurons,1);

epochs = 10; %an epoch means one full presentation % of all the characters 100 
presentationTime = 300; %Each character is presented for 300ms
timeStep = 0.2; %time step for the simulation is 0.2 ms
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
            for neuron = 1:outputNeurons
                voltageSecLayer(neuron) = neuronsTSMP(likV,weights,thresholds,Tmin,Tmax,stc,Rm);
            end
            
        end
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
end

