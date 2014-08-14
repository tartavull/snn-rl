if (debugScript)
    %Record and plot
    charMonitor.record(time,charMatrix);
    charMonitor.plot();
    meshMonitor.record(time,input);
    meshMonitor.plot();
end
if (epochIndex > 95) %To show only the last results
    debugScript = true;
end