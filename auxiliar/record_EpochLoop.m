        function obj = record_EpochLoop(obj)
            if(obj.enabled)
                obj.uniqueSpikePercentageTotal = 0;
                obj.topVsClosestNtmpTotal = 0;
                obj.topVsAllNtmpTotal = 0;
                obj.numberOfSpikesPerChar = zeros(1,4);
            endif
        endfunction 
