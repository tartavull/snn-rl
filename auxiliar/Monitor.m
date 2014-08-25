function Monitor(handle)
    %MONITOR record and plot net activity

        interactive = true;
        history;
        plotType = 'squares';
        figHandle;
        subPlot = false;
        subPlot_m;
        subPlot_n;
        subPlot_p;
        
        lineHandle;
        subPlotHandle;

        function setInteractive(obj,parameter)
            obj.interactive = parameter;
        endfunction
        function setPlotType(obj,string)
            obj.plotType = string;
        endfunction
        function setSubPlot(obj,figHandle,m,n,p)
            obj.figHandle = figHandle;
            obj.subPlot_m = m;
            obj.subPlot_n = n;
            obj.subPlot_p = p;
            obj.subPlot = true;
 
            figure(obj.figHandle);
            obj.subPlotHandle = subplot(obj.subPlot_m,obj.subPlot_n,obj.subPlot_p);
        endfunction

        function record(obj,time,data)
            obj.history{1,end+1} = time;
            obj.history{2,end} = data;
        endfunction
        
        function plot(obj)
            
            if(obj.subPlot == true && obj.subPlotHandle ~= gca)  
                subplot(obj.subPlot_m,obj.subPlot_n,obj.subPlot_p);
            endif

            switch obj.plotType
                case 'lines'
                    obj.linePlot;
                case 'lines3d'
                    obj.line3dPlot;
                case 'squares'
                    obj.squaresPlot;
                case 'char'
                    obj.charPlot;
                otherwise
                    warning('Unexpected plot type. No plot created.');
            endswitch
        endfunction
        
        function linePlot(obj)
            time = [obj.history{1,end}];
            data = [obj.history{2,end}];
            
            %Intialize the handles the first time the this is called
            if(isempty(obj.lineHandle))
                for index = 1:length(data)
                    obj.lineHandle(index) = line(nan, nan); %# Generate a blank line and return the line handle
                endfor
            endif
            
            for index = 1:length(data)
                oldTime = get(obj.lineHandle(index), 'XData');
                oldData = get(obj.lineHandle(index), 'YData');

                oldTime = [oldTime time];
                oldData = [oldData data(index)*0.8+index];
                
                timeFrame = 100;
                startIndex = max([length(oldTime) - mod(length(oldTime),timeFrame) 1]);
                set(obj.lineHandle(index), 'XData', oldTime(startIndex:end), 'YData', oldData(startIndex:end));

            endfor
            drawnow;
           
        endfunction
        
        function squaresPlot(obj)            

            if( length(obj.history) > 20)
                start = length(obj.history) - 20;
            else
                start = 1;
            endif
            
            data = [obj.history{2,start:end}];
            
            pcolor([[data zeros(size(data,1),1)] ; zeros(1,size(data,2)+1)]);
            colormap([0 0 0 ; 0 1 0]);
            axis ij;
            axis square;
            %xlim([start,length(data)]); ylim([1,16]);  % static limits
            drawnow;
        endfunction
        
        function line3dPlot(obj)
            time = [obj.history{1,end}];
            data = [obj.history{2,end}];
            
            %Intialize the handles the first time the this is called
            if(isempty(obj.lineHandle))
                colorspec = {'r'; 'g'; 'b';'c';'m';'y';'k'};
                for index = 1:length(data)
                    obj.lineHandle(index) = line(nan, nan, nan); %# Generate a blank line and return the line handle
                    set(obj.lineHandle(index),'Color',colorspec{mod(index,5)+1});
                endfor
                zlim('manual');
                zlim([0 15]);
                view(3);
                grid on;
            endif
            
            for index = 1:length(data)
                oldTime = get(obj.lineHandle(index), 'XData');
                oldPosition = get(obj.lineHandle(index), 'YData');
                oldData = get(obj.lineHandle(index), 'ZData');

                oldTime = [oldTime time];
                oldPosition = [oldPosition index];
                oldData = [oldData data(index)];
                
                timeFrame = 500;
                startIndex = max([(length(oldTime) - timeFrame) 1]);
                set(obj.lineHandle(index), 'XData', oldTime(startIndex:end), 'YData', oldPosition(startIndex:end),'ZData',oldData(startIndex:end));

            endfor
            drawnow;
        endfunction
        
        function charPlot(obj)         
           data = obj.history{2,end};
           obj.showChar(data);
        endfunction

        function handle = showChar(charArray)
            handle = imshow(charArray,'InitialMagnification',1000);
        endfunction
        function handle = plotFirings( firings )
            handle = plot(firings(:,1),firings(:,2),'.');
        endfunction

endfunction

