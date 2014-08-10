function lR = logResults(letter, epochIndex, voltagesMembraneTotal, addsDiracForChar, fid)
% Log Results: Record neuron activity results as a tsv file for analysis
%   letter = input letter 
%   epochIndex = which epoch is represented
%   voltagesMembraneTotal = total membrane voltages from neurons
%   addsDiracInEpoch = neurons that fired spikes 

fprintf(fid, '%s\t', letter);
fprintf(fid, '%d\t', epochIndex);
fprintf(fid, '%d\t', voltagesMembraneTotal);
fprintf(fid, '%d\t', addsDiracForChar);
fprintf(fid, '\n');

end
