function obj = Logger()
    %LOGGER log variables to file
    %    Record neuron activity results as a tsv file for analysis

            %Clear output and open for writing
            obj.fid = fopen('logOfResults.tsv','w');
            fclose(obj.fid);
            obj.fid = fopen('logOfResults.tsv','a');
endfunction

