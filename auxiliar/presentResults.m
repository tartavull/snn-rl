        function presentResults(obj, numberOfChars, letters)
            if(obj.enabled)
                fprintf('Avg Percent of Spikes per Char:  ');
                fprintf('%f\t', obj.uniqueSpikePercentageTotal/numberOfChars);
                fprintf('\n');
                
                fprintf('Spikes per Char:  ');
                for letterIndex = 1:size(obj.numberOfSpikesPerChar, 2)
                    fprintf('%s: %f  ', letters{letterIndex}, obj.numberOfSpikesPerChar(1, letterIndex));
                endfor
                fprintf('\n');
                
                fprintf('topVsClosestNtmp:  ');
                fprintf('%f\t', obj.topVsClosestNtmpTotal/numberOfChars);
                fprintf('\n');
                
                fprintf('topVsAllNtmp:  ');
                fprintf('%f\t', obj.topVsAllNtmpTotal/numberOfChars);
                fprintf('\n');
            endif
        endfunction
