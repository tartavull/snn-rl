debugScript = true;
if (debugScript)
    %Construct input Monitor
    meshMonitor = Monitor;
    handle = figure(1);
    meshMonitor.setSubPlot(handle,3,2,2);
    %Construct Char Monitor
    charMonitor = Monitor;
    charMonitor.setSubPlot(handle,3,2,1);
    charMonitor.setPlotType('char');
    %Construct LIK Spike Monitor
    spikeMonitor = Monitor;
    spikeMonitor.setSubPlot(handle,3,2,[3 4]);
    spikeMonitor.setPlotType('lines3d');
    %Construct adds Spike Monitor
    addsMonitor = Monitor;
    addsMonitor.setSubPlot(handle,3,2,[5 6]);
    addsMonitor.setPlotType('lines3d');
end
