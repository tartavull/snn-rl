classdef Monitor < handle
    %MONITOR record and plot net activity
    
    properties
        interactive = true;
        history;
        plotType = 'squares';
        figHandle;
        subPlot = false;
        subPlot_m;
        subPlot_n;
        subPlot_p;
        
        lHandle = 0;
    end
    
    methods
        function setInteractive(obj,parameter)
            obj.interactive = parameter;
        end
        function setPlotType(obj,string)
            obj.plotType = string;
        end
        function setSubPlot(obj,figure,m,n,p)
            obj.figHandle = figure;
            obj.subPlot_m = m;
            obj.subPlot_n = n;
            obj.subPlot_p = p;
            obj.subPlot = true;
        end
        
        
        function record(obj,time,data)
            obj.history{1,end+1} = time;
            obj.history{2,end} = data;
        end
        
        function plot(obj)
            if(obj.subPlot)  
                figure(obj.figHandle);
                subplot(obj.subPlot_m,obj.subPlot_n,obj.subPlot_p);
            end
            
            
            switch obj.plotType
                case 'lines'
                    obj.realtimePlot;
                case 'squares'
                    obj.squaresPlot;
                case 'char'
                    obj.charPlot;
                otherwise
                    warning('Unexpected plot type. No plot created.');
            end
        end
        
        function linePlot(obj)
            
            %This doesn't work yet! make something faster
            if( length([obj.history{1,:}]) > 50)
                start = length(obj.history) - 50;
            else
                start = 1;
            end
                        
            time = [obj.history{1,start:end}];
            data = [obj.history{2,start:end}];
            for line = 1:length(obj.history{2,1})
                hold on;
                plot(time,data(line,:)*0.8+line);
            end
            drawnow;
            obj.figHandle = gcf();
            
        end
        function realtimePlot(obj)
            if( obj.lHandle  == 0)
                obj.lHandle = line(nan, nan); %# Generate a blank line and return the line handle
            end
            
            time = [obj.history{1,end}];
            data = [obj.history{2,end}];
            
            X = get(obj.lHandle, 'XData');
            Y = get(obj.lHandle, 'YData');
                
            X = [X time];
            Y = [Y data];
                
            set(obj.lHandle, 'XData', X, 'YData', Y);
           
        end
        
        function squaresPlot(obj)            

            if( length(obj.history) > 20)
                start = length(obj.history) - 20;
            else
                start = 1;
            end
            
            data = [obj.history{2,start:end}];
            
            pcolor([[data zeros(size(data,1),1)] ; zeros(1,size(data,2)+1)]);
            colormap([0 0 0 ; 0 1 0]);
            axis ij;
            axis square;
            %xlim([start,length(data)]); ylim([1,16]);  % static limits
            drawnow;
            obj.figHandle = gcf();
        end
        
        function charPlot(obj)         
           data = obj.history{2,end};
           obj.showChar(data);
           obj.figHandle = gcf();

        end
    end
    
    methods(Static)
        function handle = showChar(charArray)
            handle = imshow(charArray,'InitialMagnification',1000);
            
        end
    end

    
end

