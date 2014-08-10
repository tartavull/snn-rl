function [ output_args ] = presentResults( uniqueSpikePercentageTotal, topVsClosestNtmpTotal, topVsAllNtmpTotal, numberOfChars )
%Present Results: presents a description of performance results
%   uniqueSpikePercentage = the average number of spikes produced for each 
%   Char.  It represents the amount of neurons that uniquely identify with each 
%   Char.
%   topVsClosestNtmpTotal = The difference in membrane potential between
%   the top potential and the one closest to it.  Represents how distinctly
%   the top responding neuron's potential was compared to the 2cnd highest.
%   topVsAllNtmpTotal = The difference in membrane potential between
%   the top potential and the average of the others.  Represents how distinctly
%   the top responding neuron's potential was compared to the average of the others.


fprintf('Average Percentage of Spikes per Char:');
fprintf('%d\t', uniqueSpikePercentageTotal/numberOfChars);
fprintf('\n');

fprintf('Difference in Top Highest Neuron Total Membrane Potenial (Ntmp) and Potential of Closest Ntmp:');
fprintf('%d\t', topVsClosestNtmpTotal/numberOfChars);
fprintf('\n');

fprintf('Difference in Top Highest Neuron Total Membrane Potenial (Ntmp) and Average Potential of All Other Ntmp:');
fprintf('%d\t', topVsAllNtmpTotal/numberOfChars);
fprintf('\n');

end

