function obj = CombinedMonitors()
    %COMBINEDMONITORS Use several monitors to plot the activity of a neural
    %network

        handle = figure();
        monitors = struct();
        enabled = true;

        function obj = CombinedMonitors()
            
            %Construct input Monitor
            obj.monitors.meshMonitor = Monitor;
            obj.monitors.meshMonitor.setSubPlot(obj.handle,3,2,2);
            %Construct Char Monitor
            obj.monitors.charMonitor = Monitor;
            obj.monitors.charMonitor.setSubPlot(obj.handle,3,2,1);
            obj.monitors.charMonitor.setPlotType('char');
            %Construct LIK Spike Monitor
            obj.monitors.spikeMonitor = Monitor;
            obj.monitors.spikeMonitor.setSubPlot(obj.handle,3,2,[3 4]);
            obj.monitors.spikeMonitor.setPlotType('lines3d');
            %Construct adds Spike Monitor
            obj.monitors.addsMonitor = Monitor;
            obj.monitors.addsMonitor.setSubPlot(obj.handle,3,2,[5 6]);
            obj.monitors.addsMonitor.setPlotType('lines3d');
            
        endfunction
        
        function record_EpochLoop(obj)
            if(obj.enabled)
            endif
        endfunction
        
        function record_DictionaryLoop(obj,time,charMatrix,input)
            if(obj.enabled)
                obj.monitors.charMonitor.record(time,charMatrix);
                obj.monitors.meshMonitor.record(time,input);
                
                obj.monitors.charMonitor.plot();
                obj.monitors.meshMonitor.plot();
            endif
        endfunction
        
        function record_StepLoop(obj,time,likV,voltagesMembrane)
            if(obj.enabled)
                obj.monitors.spikeMonitor.record(time,likV);
                obj.monitors.addsMonitor.record(time,voltagesMembrane);
                
                obj.monitors.spikeMonitor.plot();
                obj.monitors.addsMonitor.plot();
            endif
        endfunction
    
endfunction
