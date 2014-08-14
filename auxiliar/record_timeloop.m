
if (debugScript)
    spikeMonitor.record(time,likV);
    spikeMonitor.plot();
    addsMonitor.record(time,voltagesMembrane);
    addsMonitor.plot();
end