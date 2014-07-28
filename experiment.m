close all;
clear all;
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

epochs = 5; %an epoch means one full presentation % of all the characters
presentationTime = 300; %Each character is presented for 300ms
timeStep = 0.2; %time step for the simulation is 0.2 ms

%Construct Spike Monitor
monitor = Monitor;
handle = figure(1);
monitor.setSubPlot(handle,1,2,2);
plotFrecuency = 10; %The plot are updated every 10ms


index=0;
for time = 0:timeStep:presentationTime*length(Dictionary)*epochs
    
    if(mod(time,presentationTime)==0)
        %Each character is presented one at a time sequentially during
        %the training process
        charCounter= mod(time/presentationTime,length(Dictionary))+1;
        
        %Plot new char
         charMatrix = Dictionary{charCounter,2};
         subplot(1,2,1);
         hold on;
         imshow(charMatrix,'InitialMagnification',1000);
         drawnow;
        

        %Update input
        input = reshape(charMatrix,[],1);
        
        %Print epoch number
        if(mod(time,presentationTime*length(Dictionary))==0)
            epochCounter=time/(presentationTime*length(Dictionary));
            fprintf('Epoch %d of %d \n',epochCounter,epochs);
        end
    end
    
    if(mod(time,plotFrecuency)==0)
        monitor.record(time,input);
    end
    
end