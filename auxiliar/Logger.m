classdef Logger < handle
    %LOGGER log variables to file
    %    Record neuron activity results as a tsv file for analysis
    
    properties
        enabled = true;
        fid
    end
    
    methods
        function obj = Logger(enabled)
            obj.enabled = enabled;
            if(enabled)
                %Clear output and open for writing
                obj.fid = fopen('logOfResults.tsv','w');
                fclose(obj.fid);
                obj.fid = fopen('logOfResults.tsv','a');
            end
        end
        function delete(obj)
            if(enabled)
                fclose(obj.fid);
            end
        end
        function record_DictionaryLoop(obj,letter,epochIndex,voltagesMembraneTotal,addsDiracForChar)
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
            end
        end
    end
    
end

