function [ percentageOfUniqueSpikes, topVsClosestNtmp, topVsAllNtmp ] = scoreResults( letter, epochIndex, voltagesMembraneTotal, addsDiracForChar )
%Score results: Calculates values used in scoring

percentageOfUniqueSpikes = size(find(addsDiracForChar == 1), 2)/size(addsDiracForChar, 2);

sortedVoltageMembraneTotals = sort(voltagesMembraneTotal);

topVsClosestNtmp = sortedVoltageMembraneTotals(end)-sortedVoltageMembraneTotals(end-1);

topVsAllNtmp = sortedVoltageMembraneTotals(end)-(sum(sortedVoltageMembraneTotals(1:end-1))/(size(sortedVoltageMembraneTotals, 2) - 1));

end

