function [ output_args ] = presentResults( uniqueSpikePercentageTotal, numberOfSpikesPerChar, Dictionary, topVsClosestNtmpTotal, topVsAllNtmpTotal, numberOfChars )
%Present Results: presents a description of performance results
%   uniqueSpikePercentage = the average number of spikes produced for each 
%   Char.  It represents the amount of neurons that uniquely identify with each 
%   Char.
%   numberOfSpikesPerChar = Spikes fired for each Char
%   topVsClosestNtmpTotal = The difference in membrane potential between
%   the top potential and the one closest to it.  Represents how distinctly
%   the top responding neuron's potential was compared to the 2cnd highest.
%   topVsAllNtmpTotal = The difference in membrane potential between
%   the top potential and the average of the others.  Represents how distinctly
%   the top responding neuron's potential was compared to the average of the others.


fprintf('Avg Percent of Spikes per Char:  ');
fprintf('%d\t', uniqueSpikePercentageTotal/numberOfChars);
fprintf('\n');

fprintf('Spikes per Char:  ');
for letterIndex = 1:size(numberOfSpikesPerChar, 2)
    fprintf('%s: %d  ', Dictionary{letterIndex}, numberOfSpikesPerChar(1, letterIndex));
end
fprintf('\n');

fprintf('topVsClosestNtmp:  ');
fprintf('%d\t', topVsClosestNtmpTotal/numberOfChars);
fprintf('\n');

fprintf('topVsAllNtmp:  ');
fprintf('%d\t', topVsAllNtmpTotal/numberOfChars);
fprintf('\n');

end

