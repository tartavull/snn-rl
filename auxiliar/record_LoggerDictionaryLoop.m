        function obj = record_LoggerDictionaryLoop(obj,letter,epochIndex,voltagesMembraneTotal,addsDiracForChar)
            if(obj.enabled)
                %   letter = input letter
                %   epochIndex = which epoch is represented
                %   voltagesMembraneTotal = total membrane voltages from neurons
                %   addsDiracInEpoch = neurons that fired spikes
                
                fprintf(obj.fid, '%s\t', letter);
                fprintf(obj.fid, '%d\t', epochIndex);
                fprintf(obj.fid, '%d\t', voltagesMembraneTotal);
                fprintf(obj.fid, '%d\t', addsDiracForChar);
                fprintf(obj.fid, '\n');
            endif
        endfunction
