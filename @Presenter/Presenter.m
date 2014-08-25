function presenter = Presenter()
    %PRESENTER  presents a description of performance results
    %   uniqueSpikePercentage = the average number of spikes produced for each
    %   Char.  It represents the amount of neurons that uniquely identify
    %   with each Char.
    %   numberOfSpikesPerChar = Spikes fired for each Char
    %   topVsClosestNtmpTotal = The difference in membrane potential between
    %   the top potential and the one closest to it.  Represents how distinctly
    %   the top responding neuron's potential was compared to the 2cnd highest.
    %   topVsAllNtmpTotal = The difference in membrane potential between
    %   the top potential and the average of the others.  Represents how distinctly
    %   the top responding neuron's potential was compared to the average of the others.
    

        presenter.enabled = true;
        
        presenter.uniqueSpikePercentageTotal = 0;
        presenter.topVsClosestNtmpTotal = 0;
        presenter.topVsAllNtmpTotal = 0;
        presenter.numberOfSpikesPerChar = zeros(1,4);

endfunction
        
