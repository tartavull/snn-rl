        function obj = record_PresenterDictionaryLoop(obj, dictionaryIndex , voltagesMembraneTotal, addsDiracForChar);
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
                
            endif
        endfunction
