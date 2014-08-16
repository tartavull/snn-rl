close all;
clear all;
addpath('adds','datasets','auxiliar');

%Initialize the net
architecture;

%Initialize monitors for plotting
combinedMonitors = CombinedMonitors();
%combinedMonitors.enabled = false;

%Save variables of interest to disk
logger = Logger();
logger.enabled = false;

%Print info about net performance
presenter = Presenter();
presenter.enabled = false;


time = 0;
for epochIndex = 1:epochs
    
    presenter.record_EpochLoop();
    
    for dictionaryIndex = 1:length(Dictionary)
        %Each character is presented one at a time sequentially during
        %the training process
        charCounter= mod(dictionaryIndex,length(Dictionary))+1;
        charMatrix = Dictionary{charCounter,2};
        input = reshape(charMatrix,[],1);
        
        %Used to print useful information
        voltagesMembraneTotal = 0;
        addsDiracForChar=zeros(1,4);
        
        %Used to plot current character
        combinedMonitors.record_DictionaryLoop(time,charMatrix,input);
        
        %Main code
        for stepIndex = 1:presentationTime*timeStep
            time = time + timeStep;
            
            %Integrate and fire layer
            likI = input * stimulusIntensity;
            likFired=find(likV > firingThreshold);    % indices of spikes
            likV = likV + timeStep  .* 1/capacitance * (likI - likV./ resistance);
            likV(likFired) = restPotential; %Set to resting potential
            
            likDirac = zeros(15,1);
            likDirac(likFired) = 1;
            likFirings(:,stepIndex) = likDirac; %Something isn't working here, almost all zeros check Monitor.plotFirings
            likLastTimeFired = [likLastTimeFired likLastTimeFired(:,end)];
            likLastTimeFired(likFired,end) = time;
            
            %Second layer
            addsFired=find(voltagesMembrane > firingThreshold);
            addsDirac = zeros(1,4);
            addsDirac(addsFired) = 1;
            addsFirings(:,stepIndex) = addsDirac;
            addsLastTimeFired = [addsLastTimeFired addsLastTimeFired(:,end)];
            addsLastTimeFired(addsFired,end) = time;
            
            %For logger and presenter
            addsDiracForChar(addsFired) = 1;
            
            tauDendritic = timeConstant(tauMax, tauMin , weightsDendritic);
            resistenceDendritic = resistanceComputation(tauDendritic, firingThreshold , resistenceMembrane, tauMembrane);
            
            currentDendritic = dendriticPostSynapticCurrent(timeStep, currentDendritic , resistenceDendritic, weightsDendritic, tauDendritic, likDirac);
            currentSomatic = currentSomatic + timeStep * ((-currentSomatic + sum(weightsSomatic .* [addsDirac; addsDirac; addsDirac; addsDirac],1))./tauSomatic);
            
            voltagesMembrane =  voltagesMembrane + timeStep * ((-voltagesMembrane + resistenceMembrane .* ( sum(currentDendritic,1) + currentSomatic))/tauMembrane);
            voltagesMembrane(addsFired) = restPotential;
            
            %updateWeights
            for addsNeuron = 1:length(Dictionary)
                deltaSynapticSpike = addsLastTimeFired(addsNeuron,end)*ones(size(likLastTimeFired(end))) -  likLastTimeFired(:,end);
                deltaSynapticWeight = deltaWeight(deltaSynapticSpike);
                weightsDendritic(:,addsNeuron) = newWeight(deltaSynapticWeight,weightsDendritic(:,addsNeuron), weightMinExcitatory, weightMaxExcitatory,learningRate);
                
                somaticSynapseWeigthIndexes = setdiff(1:length(Dictionary),addsNeuron);
                deltaSomaticSpike = addsLastTimeFired(addsNeuron,end)*ones(length(Dictionary)-1,1) -  addsLastTimeFired(somaticSynapseWeigthIndexes,end);
                deltaSomaticWeight = deltaWeight(deltaSomaticSpike);
                weightsSomatic(somaticSynapseWeigthIndexes,addsNeuron) = newWeight(deltaSomaticWeight, weightsSomatic(somaticSynapseWeigthIndexes,addsNeuron), weightMinInhibitory, weightMaxInhibitory,learningRate);
            end
            
            %Used to plot membrane potentials
            combinedMonitors.record_StepLoop(time,likV,voltagesMembrane);
            
            %Used to print useful information
            voltagesMembraneTotal = voltagesMembraneTotal + voltagesMembrane;
        end
        
        
        %Used to print useful information
        logger.record_DictionaryLoop(Dictionary{charCounter,1}, epochIndex, voltagesMembraneTotal, addsDiracForChar);
        presenter.record_DictionaryLoop(dictionaryIndex, voltagesMembraneTotal, addsDiracForChar);
        
    end
    %Print epoch number
    fprintf('Epoch %d of %d \n',epochIndex,epochs);
    
    %Print selected results
    if (epochIndex == 1 || epochIndex == 2 || epochIndex == 3 || epochIndex == 25 || epochIndex == 50 || epochIndex == 75 || epochIndex == 100)
        presenter.presentResults(length(Dictionary),Dictionary(:,1));
    end
    
end
