close all;
clear all;
makeDictionary;
Dictionary = Dictionary(1:4,:); % Using A,B,C,D

epochs = 20; %Total number of epochs, (an epoch means one full presentation
             %of all the four characters) 

presentationTime = 300; %Each character is presented for 300ms

timeStep = 0.2; %time step for the simulation is 0.2 ms

graphResolution = 20; %number of states displayed in the inputHistoy plot

input = zeros(15,1);

inputHistory = zeros(15,presentationTime*length(Dictionary)*epochs/timeStep);
graphInputHistory = zeros(15,length(Dictionary)*epochs);

timeHistory = zeros(1,presentationTime*length(Dictionary)*epochs/timeStep);
    
index = 0;
time = 0;
state = 0;
for epochCounter = 1:epochs
    
    fprintf('Epoch %d of %d \n',epochCounter,epochs);        %Print epoch number
    tic;
    
    for charCounter = 1:length(Dictionary)
        %Each character is presented one at a time sequentially during
        %the training process
        
        state = state + 1;
        
        %Plot new char
        charMatrix = Dictionary{charCounter,2};
        input = reshape(charMatrix,[],1);
        subplot(1,2,2);
        imshow(charMatrix,'InitialMagnification',1000);
        
        %Save input in each change of state, just for graphing
        graphInputHistory(:,state) = input;
            
        subplot(1,2,1);
        pcolor([[graphInputHistory zeros(size(graphInputHistory,1),1)] ; zeros(1,size(graphInputHistory,2) + 1) ] )
        %plotting only the last "graphResolution" states
        colormap([0 0 0 ; 0 1 0]);
        axis ij;
        axis square;
        xlim([(state-graphResolution-1)*(state>graphResolution)+1,state+1]); ylim([1,16]);  % static limits
        drawnow;
        
        %% Faster MATRIX version, suitable for calculation but incompatible with a "live" changing graph
%         a = (epochCounter-1)*length(Dictionary)*presentationTime/timeStep + (charCounter-1)*presentationTime/timeStep + 1;
%         b = (epochCounter-1)*length(Dictionary)*presentationTime/timeStep + (charCounter)*presentationTime/timeStep;
%         
%         timeHistory(a:b) = a*timeStep : timeStep : b*timeStep;
%         inputHistory(:,a:b) = input * ones(presentationTime/timeStep,1)';

        %% Slow FOR version, suitable for debugging/displaying Voltages and Currents graphs over time
       for msCounter = 1:presentationTime/timeStep
           time = time + timeStep;
           index = index + 1;

           %Save input in each ms
           timeHistory(index) = time;
           inputHistory(:,index) = input;
       end

    end

    toc;
end