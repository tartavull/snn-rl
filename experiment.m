close all;
clear all;
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

epochs = 1; %Total number of epochs, (an epoch means one full presentation
% of all the four characters) 

presentationTime = 300; %Each character is presented for 300ms

timeStep = 0.2; %time step for the simulation is 0.2 ms

input = zeros(15,1);

plotFrecuency = 10; %The plot are updated every 10ms
inputHistory = zeros(15,presentationTime*length(Dictionary)*epochs/plotFrecuency);
timeHistory = zeros(1,presentationTime*length(Dictionary)*epochs/plotFrecuency);;

index=0;
for time = 0:timeStep:presentationTime*length(Dictionary)*epochs
    
    if(mod(time,presentationTime)==0)
        %Each character is presented one at a time sequentially during
        %the training process
        charCounter= mod(time/presentationTime,length(Dictionary))+1;
        
        %Plot new char
        charMatrix = Dictionary{charCounter,2};
        subplot(2,2,3);
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
        index=index+1;

        %Plot input
        timeHistory(index)=time;
        inputHistory(:,index) = input;
        subplot(2,1,1);
        for line = 1:length(input)
            hold on;
            plot(timeHistory(1:index),inputHistory(line,1:index)*line);
        end
        drawnow;
    end
    
end

