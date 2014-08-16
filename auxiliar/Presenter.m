classdef Presenter < handle
    %PRESENTER  presents a description of performance results
    %   uniqueSpikePercentage = the average number of spikes produced for each
    %   Char.  It represents the amount of neurons that uniquely identify
    %   with each Char.
    %   numberOfSpikesPerChar = Spikes fired for each Char
    %   topVsClosestNtmpTotal = The difference in membrane potential between
    %   the top potential and the one closest to it.  Represents how distinctly
    %   the top responding neuron's potential was compared to the 2cnd highest.
    %   topVsAllNtmpTotal = The difference in membrane potential between
    %   the top potential and the average of the others.  Represents how distinctly
    %   the top responding neuron's potential was compared to the average of the others.
    
    properties
        enabled = true;
        
        uniqueSpikePercentageTotal;
        topVsClosestNtmpTotal;
        topVsAllNtmpTotal;
        numberOfSpikesPerChar;
        
    end
    
    methods
        function record_EpochLoop(obj)
            if(obj.enabled)
                obj.uniqueSpikePercentageTotal = 0;
                obj.topVsClosestNtmpTotal = 0;
                obj.topVsAllNtmpTotal = 0;
                obj.numberOfSpikesPerChar = zeros(1,4);
            end
        end
        function record_DictionaryLoop(obj, dictionaryIndex , voltagesMembraneTotal, addsDiracForChar);
            if(obj.enabled)
                percentageOfUniqueSpikes = size(find(addsDiracForChar == 1), 2)/size(addsDiracForChar, 2);
                
                sortedVoltageMembraneTotals = sort(voltagesMembraneTotal);
                
                topVsClosestNtmp = sortedVoltageMembraneTotals(end)-sortedVoltageMembraneTotals(end-1);
                
                topVsAllNtmp = sortedVoltageMembraneTotals(end)-(sum(sortedVoltageMembraneTotals(1:end-1))/(size(sortedVoltageMembraneTotals, 2) - 1));
                
                obj.uniqueSpikePercentageTotal = obj.uniqueSpikePercentageTotal + percentageOfUniqueSpikes;
                obj.topVsClosestNtmpTotal = obj.topVsClosestNtmpTotal + topVsClosestNtmp;
                obj.topVsAllNtmpTotal = obj.topVsAllNtmpTotal + topVsAllNtmp;
                
                % Record spikes for each Char
                numberOfSpikesForChar = size(find(addsDiracForChar == 1), 2);
                obj.numberOfSpikesPerChar(1, dictionaryIndex) = numberOfSpikesForChar;
                
            end
        end
        function record_StepLoop(obj)
            if(obj.enabled)
            end
        end
        function presentResults(obj, numberOfChars, letters)
            if(obj.enabled)
                fprintf('Avg Percent of Spikes per Char:  ');
                fprintf('%d\t', obj.uniqueSpikePercentageTotal/numberOfChars);
                fprintf('\n');
                
                fprintf('Spikes per Char:  ');
                for letterIndex = 1:size(obj.numberOfSpikesPerChar, 2)
                    fprintf('%s: %d  ', letters{letterIndex}, obj.numberOfSpikesPerChar(1, letterIndex));
                end
                fprintf('\n');
                
                fprintf('topVsClosestNtmp:  ');
                fprintf('%d\t', obj.topVsClosestNtmpTotal/numberOfChars);
                fprintf('\n');
                
                fprintf('topVsAllNtmp:  ');
                fprintf('%d\t', obj.topVsAllNtmpTotal/numberOfChars);
                fprintf('\n');
            end
        end
        
        
    end
end
