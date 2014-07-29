close all;
clear all;
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

%Construct Spike Monitor
monitor = Monitor;
handle = figure(1);
monitor.setSubPlot(handle,1,2,2);
%Construct Char Monitor
charMonitor = Monitor;
charMonitor.setSubPlot(handle,1,2,1);
charMonitor.setPlotType('char');

epochs = 50; %an epoch means one full presentation % of all the characters
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
            

            
        end
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
end

