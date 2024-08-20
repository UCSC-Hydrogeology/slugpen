%%% =======================================================================
%%  Purpose:
%     This function adjusts the layout of plots in SlugPen. This includes
%     turning on and off depth, tilt, and acceleration plots. Temperature
%     plot cannot be turned off.
%%  Last edit
%     07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function plotLayout(H,...
                    DepthPlot, ...
                    AccelPlot, ...
                    TiltPlot, ...
                    grid_PlotWindow)
            
       %% Turn all values of plot checkboxes to strings
            dp = num2str(DepthPlot);
            ap = num2str(AccelPlot);
            tp = num2str(TiltPlot);
           
      %% Concatenate strings of values with the letter x so that this
      %% string can be used in the grid layout function
            dpx = strcat(dp, 'x');
            apx = strcat(ap, 'x');
            tpx = strcat(tp, 'x');

            grid_PlotWindow.RowHeight = {'0.2x','0.27x','5x', dpx, apx, tpx, '0.2x', '0.2x'};

      %% Assign all children of each axes to variable
            DC = H.Axes.Depth.Children;
            AC = H.Axes.Accel.Children;
            TC = H.Axes.Tilt.Children;

      %% Adjust grid with plots based on which plots are turned on and off
            
            % If depth, acceleration, and tilt plots are all turned off
            if DepthPlot == 0 && AccelPlot == 0 && TiltPlot == 0
                H.Axes.Depth.Visible = 'off';
                C = DC;
                axesChildrenOFF(C)

                H.Axes.Accel.Visible = 'off';
                C = AC;
                axesChildrenOFF(C)

                H.Axes.Tilt.Visible = 'off';
                C = TC;
                axesChildrenOFF(C)

            % If depth and tilt plots are off but acceleration is on
            elseif DepthPlot == 0 && AccelPlot == 1 && TiltPlot == 0
                H.Axes.Depth.Visible = 'off';
                C = DC;
                axesChildrenOFF(C)

                H.Axes.Accel.Visible = 'on';
                C = AC;
                axesChildrenON(C)

                H.Axes.Tilt.Visible = 'off';
                C = TC;
                axesChildrenOFF(C)

            % If depth and acceleration plots are off but tilt is on
            elseif DepthPlot == 0 && AccelPlot == 0 && TiltPlot == 1
                H.Axes.Depth.Visible = 'off';
                C = DC;
                axesChildrenOFF(C)

                H.Axes.Accel.Visible = 'off';
                C = AC;
                axesChildrenOFF(C)

                H.Axes.Tilt.Visible = 'on';
                C = TC;
                axesChildrenON(C)

            % If tilt and acceleration plots are off but depth is on
            elseif DepthPlot == 1 && AccelPlot == 0 && TiltPlot == 0
                H.Axes.Depth.Visible = 'on';
                C = DC;
                axesChildrenON(C)

                H.Axes.Accel.Visible = 'off';
                C = AC;
                axesChildrenOFF(C)

                H.Axes.Tilt.Visible = 'off';
                C = TC;
                axesChildrenOFF(C)

            % If tilt and acceleration plots are on but tilt is off
            elseif DepthPlot == 1 && AccelPlot == 1 && TiltPlot == 0
                H.Axes.Depth.Visible = 'on';
                C = DC;
                axesChildrenON(C)

                H.Axes.Accel.Visible = 'on';
                C = AC;
                axesChildrenON(C)

                H.Axes.Tilt.Visible = 'off';
                C = TC;
                axesChildrenOFF(C)

            % If tilt and depth plots are on but acceleration is off
            elseif DepthPlot == 1 && AccelPlot == 0 && TiltPlot == 1
                H.Axes.Depth.Visible = 'on';
                C = DC;
                axesChildrenON(C)

                H.Axes.Accel.Visible = 'off';
                C = AC;
                axesChildrenOFF(C)

                H.Axes.Tilt.Visible = 'on';
                C = TC;
                axesChildrenON(C)

            % If tilt and acceleration plots are on but depth is off
            elseif DepthPlot == 0 && AccelPlot == 1 && TiltPlot == 1
                H.Axes.Depth.Visible = 'off';
                C = DC;
                axesChildrenOFF(C)

                H.Axes.Accel.Visible = 'on';
                C = AC;
                axesChildrenON(C)

                H.Axes.Tilt.Visible = 'on';
                C = TC;
                axesChildrenON(C)

            % If depth, acceleration, and tilt plots are all turned on
            elseif DepthPlot == 1 && AccelPlot == 1 && TiltPlot == 1
                H.Axes.Depth.Visible = 'on';
                C = DC;
                axesChildrenON(C)

                H.Axes.Accel.Visible = 'on';
                C = AC;
                axesChildrenON(C)

                H.Axes.Tilt.Visible = 'on';
                C = TC;
                axesChildrenON(C)
            end


