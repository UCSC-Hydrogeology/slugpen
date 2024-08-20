%%% =======================================================================
%%  Purpose: 
%       This function aligns the inner position left edge and width 
%       (x axis) of all 4 plots of SlugPen.
%%  Last edit:
%       07/18/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function [Checkboxes,TempPlots, BadSensors]...
            = PlotData(figure_Main, ...
             H, ...
             DATA, ...
             panel_Thermistors)

%% Loading dialog box
           loading = uiprogressdlg(figure_Main, 'Message', 'Plotting raw data...');
           loading.Value = .6;
       
%% Initialize 

            % Determine Number of Thermistors
            % ----------------------------------
            NoTherm = H.Fileinfo.No_Thermistors.Value;

            % Extract Decimated Data for plotting
            % -----------------------------------
            T_dec          = DATA.Tdec;
            T_decRaw       = DATA.TdecRaw;
            Depth_dec      = DATA.Depth_dec;
            G_dec          = DATA.G_dec;
            Tilt_dec       = DATA.Tilt_dec;
            Time_dec       = DATA.Time_dec;

            % Extract Clean Data for plotting
            % ---------------------------------
            T_cln          = DATA.Tcln;

            % Extract Raw Data for plotting
            % -----------------------------
            T      = DATA.Traw;
            Time   = DATA.Time;


            % Round data to 3 digits
            % ----------------------
             round(T_cln,3);
             round(T_dec,3);
             round(T,3);
             round(T_decRaw,3);

            % Enable buttons and fill in text
            % --------------------------------
            H.Selections.Select.Enable='on';
            H.Exe_Controls.Delete.Enable='on';

            % Generate colormap
            % -------------------
            colormap(H.Axes.Raw, 'default');
            CMap= colormap(H.Axes.Raw, 'turbo');
            CMap = flipud(CMap);
            CMap = interp1(1:256,CMap,1:256/(NoTherm-1):256); % 'jet' colormap for Matlab version 2019a and later has size 256x3 RBG matrix
            colormap(H.Axes.Raw, CMap);
            
            % Initialize variables
            % --------------------
            BadSensors=[];
            %h_axDec = repmat([], 1, NoTherm);
            %h_axCln = repmat([], 1, NoTherm);
            %h_axRaw = repmat([], 1, NoTherm);
            %h_axRawDec = repmat([], 1, NoTherm);



%% Plot FILTERED DECIMATED (median filtered) temperature and ohm data as o's

            % Vector to store plot handles for each thermistor 
            % ------------------------------------------------
            leg = zeros(1,NoTherm);
            i=1;
            while i<=NoTherm-1
                % Plot by time
                h_axDec(i) = plot(H.Axes.Raw, Time_dec,T_dec(i,:),'-o','markersize',3,...
                    'color',CMap(i,:),'markerfacecolor',CMap(i,:));    
                hold(H.Axes.Raw, 'on');
                leg(i) = i;
                i=i+1;
            end

            % Plot DECIMATED bottom water
            % ---------------------------
            h_axDec(i) = plot(H.Axes.Raw, Time_dec,T_dec(NoTherm,:),'k-x','markersize',3);

            % Save plot handles
            % -----------------
            H.Plot_Controls.DecimatedPlot.UserData = h_axDec;

            % Enable plot control toggle
            % ---------------------------
            H.Plot_Controls.DecimatedPlot.Enable='on';
            
            % Define decimated plot
            % ----------------------
            H.Plot_Controls.DecimatedPlot.Value = H.PlotCheckboxes.Filtereddec_check.Value; 

            % Labels, date ticks, and title
            % ---------------------------------
            if T_dec(5,10)>100
                ylabel(H.Axes.Raw, 'Resistance (Ohm)');
            end
            
            % Keep visible or hide
            % --------------------
            if H.Plot_Controls.DecimatedPlot.Value == 0
                for i=1:NoTherm
                  set(H.Plot_Controls.DecimatedPlot.UserData(i),'Visible','off');
                end
            end


%% Plot FILTERED Data
            % Plot CLEAN (median filtered) data as x's
            % Vector to store plot handles for each thermistor 

            i=1;
            while i<=NoTherm-1
                h_axCln(i) = plot(H.Axes.Raw, Time,T_cln(i,:),'o','markersize',3,...
                    'color',CMap(i,:),'markerfacecolor',CMap(i,:));
                hold(H.Axes.Raw, "on");
                i=i+1;
            end
            
            % Plot CLEAN bottom water
            % --------------------
            h_axCln(i) = plot(H.Axes.Raw, Time,T_cln(NoTherm,:),'kv','markersize',3,'markerfacecolor','b');

            % Turn bottom water off;
            % --------------------
            set(h_axCln(i),'Visible','off')


            % Save plot handles
            % --------------------
            H.Plot_Controls.CleanPlot.UserData = h_axCln;

            % Enable plot control toggle
            % --------------------------
            H.Plot_Controls.CleanPlot.Enable='on';

            % Keep visible or hide
            % --------------------
            if H.Plot_Controls.CleanPlot.Value == 0
                for i=1:NoTherm
                    set(H.Plot_Controls.CleanPlot.UserData(i),'Visible','off');
                end
    
            end

            % Hide clean data
            H.PlotCheckboxes.Filtered_check.Value = 0;

      %% RAW DATA PLOT

            % Plot RAW data as x's
            i=1;
            while i<=NoTherm-1
                h_axRaw(i) = plot(H.Axes.Raw, Time,T(i,:),'*','markersize',4,'color',CMap(i,:));
                hold(H.Axes.Raw, "on");
                i=i+1;
            end
            % Plot RAW Bottom Water
            h_axRaw(i) = plot(H.Axes.Raw, Time,T(NoTherm,:),'kv','markersize',4);


            % Save plot handles
            H.Plot_Controls.RawPlot.UserData = h_axRaw;

            % Enable Plot Control Toggle
            H.Plot_Controls.RawPlot.Enable='on';
            H.Plot_Controls.RawPlot.Value=1;

            
            H.Axes.Raw.PositionConstraint = 'innerposition'; % constrains axes box to inner position

      %% RAW DECIMATED DATA PLOT
             
            leg = zeros(1,NoTherm);
            i=1;
            while i<=NoTherm-1
                % Plot by time
                h_axRawDec(i) = plot(H.Axes.Raw, Time_dec,T_decRaw(i,:),'-o','markersize',3,...
                    'color',CMap(i,:),'markerfacecolor',CMap(i,:));    
                hold(H.Axes.Raw, 'on');
                leg(i) = i;
                i=i+1;
            end

            % Plot DECIMATED Bottom Water
            h_axRawDec(i) = plot(H.Axes.Raw, Time_dec,T_decRaw(NoTherm,:),'k-x','markersize',4);

            % Save plot handles
            H.Plot_Controls.RawDecimatedPlot.UserData = h_axRawDec;

            % Enable Plot Control Toggle
            H.Plot_Controls.RawDecimatedPlot.Enable='on';
            
            % Define decimated plot
            H.Plot_Controls.RawDecimatedPlot.Value = H.PlotCheckboxes.Raw_check.Value; 

      %% DEPTH - DECIMATED
            % Plot depth data 
            h_axDep=plot(H.Axes.Depth, Time_dec,-Depth_dec,'b');

            % Add the axes data to the Structure
            H.Axes.Depth.UserData = h_axDep;

            H.Axes.Depth.PositionConstraint = 'innerposition'; % constrains axes box to inner position

     %% ACCELERATION - DECIMATED
            % Plot Acceleration data 
            h_axAcc=plot(H.Axes.Accel, Time_dec,G_dec,'k', 'Visible','off');
            h_axAcc.Visible="off";

            % Add the axes data to the Structure
            H.Axes.Accel.UserData = h_axAcc;

            H.Axes.Accel.PositionConstraint = 'innerposition'; % constrains axes box to inner position


     %% TILT - DECIMATED
            % Plot Tilt
            plot(H.Axes.Tilt, Time_dec,Tilt_dec,'r'); hold(H.Axes.Tilt, "on")

            H.Axes.Tilt.PositionConstraint = 'innerposition'; % constrains axes box to inner position

     %% Link Axes 
            tempAx = H.Axes.Raw;
            accelAx = H.Axes.Accel;
            depthAx = H.Axes.Depth;
            tiltAx = H.Axes.Tilt;
            ax=[tempAx depthAx accelAx tiltAx];
            linkaxes(ax,'x');
            [ax.XLimMode] = deal('auto');
    
     %% Axes ticks 

            % X axis (time)
            t_start = Time_dec(1);
            t_end = Time_dec(end);
            tempAx.XLim = [t_start,t_end];
            H.Axes.Raw.XLim=tempAx.XLim;
            H.Axes.Depth.XLim=tempAx.XLim;
            H.Axes.Accel.XLim=tempAx.XLim;
            H.Axes.Tilt.XLim=tempAx.XLim;

            tempAx.XAxis.TickLabelFormat='HH:mm:ss';
            accelAx.XAxis.TickLabelFormat='HH:mm:ss';
            depthAx.XAxis.TickLabelFormat='HH:mm:ss';
            tiltAx.XAxis.TickLabelFormat='HH:mm:ss';

            tempAx.XTickMode='auto';
            accelAx.XTickMode='auto';
            depthAx.XTickMode='auto';
            tiltAx.XTickMode='auto';

            tempAx.XAxis.SecondaryLabel.Visible='off';
            accelAx.XAxis.SecondaryLabel.Visible='off';
            depthAx.XAxis.SecondaryLabel.Visible='off';
            tiltAx.XAxis.SecondaryLabel.Visible='off';

            tempAx.XMinorTick='on';
            accelAx.XMinorTick='on';
            depthAx.XMinorTick='on';
            tiltAx.XMinorTick='on';

            % Y axis

            tempAx.YLim = [min(min(T))-1,max(max(T)+1)];

            H.Axis_Limits.AxLims_Raw_Y=tempAx.YLim;

    

    %% Update loading box
            loading.Message = 'Finishing...';
            loading.Value = .9;

    %% Plot sensor checkboxes
    %% ----------------------

            NumberOfSensors = NoTherm;
            
                % Grid layout for sensor panel
                % -----------------------------
                    
                    for i=NumberOfSensors:-1:1
                        ThermistorGrid = uigridlayout(panel_Thermistors, 'BackgroundColor', ...
                        [0.94 0.94 0.94], 'ColumnWidth', {'100x'}, 'RowHeight', cellstr(repmat('1x', NumberOfSensors, 1))');
                    end
            
            
                % Create panels for each sensor
                % ------------------------------
                    for i=NumberOfSensors:-1:1
                        SensorPanels{i} = uipanel(ThermistorGrid, 'BorderType', ...
                            'none', 'TitlePosition', 'centertop', ...
                            'Title', [' T' num2str(i)], 'FontWeight', 'bold', ...
                            'BackgroundColor', [0.94 0.94 0.94]);
                    end
            
            
                % Grid layout for each sensor
                % -----------------------------
                    for i=NumberOfSensors-1:-1:1
                        SensorGrids{i} = uigridlayout(SensorPanels{i}, 'ColumnWidth', ...
                            {'100x'}, 'RowHeight', {'100x'}, 'Padding',[0 0 0 0], ...
                            'BackgroundColor', CMap(i,:));
                    end
                    i = max(NumberOfSensors);
                    SensorGrids{i} = uigridlayout(SensorPanels{i}, 'ColumnWidth', ...
                            {'100x'}, 'RowHeight', {'100x'}, 'Padding',[0 0 0 0], ...
                            'BackgroundColor', 'k');
                    SensorPanels{i}.Title = 'BW';
            
                % Checkbox for each sensor
                % -----------------------------
            
                  for i=NumberOfSensors:-1:1
                          Checkboxes{i} = uicheckbox(SensorGrids{i}, 'Text', '', ...
                              'Value',true, 'Tooltip', ['Check to display data from sensor T' num2str(i) ' on plots'], ...
                              'tag',['c_b_' num2str(i)], 'ValueChangedFcn', {@cbSensValueChange, h_axDec(i), h_axCln(i), h_axRaw(i), h_axRawDec(i) ...
                              H.PlotCheckboxes.Filtered_check.Value, H.PlotCheckboxes.Filtereddec_check.Value, H.PlotCheckboxes.Raw_check.Value, ...
                              H.PlotCheckboxes.Rawdec_check.Value});
                  end
     
     %% Save plots in a structure
        TempPlots = struct('FilteredDecimated',h_axDec, ...
            'Filtered', h_axCln, 'RawDecimated', h_axRawDec,...
            'Raw',h_axRaw);


     %% Hide plots and checkboxes for bad (NaN) data
        
        % Format compatibility
        T(isnan(T))=-999;

        if any(all(T==-999,2))
            hidedata = uiconfirm(figure_Main, ['One of more sensor recorded ' ...
                'no data (T = -999) in this .raw file. Would you like to ' ...
                'hide this data form the plot?'], 'NaN Data Recorded', ...
                'Options',{'Hide NaN (T = -999) data', 'Plot NaN (T = -999) data'});
            switch hidedata
                case 'Hide NaN (T = -999) data'
                    BadSensors = find(all(T==-999,2));
                    for i=1:length(BadSensors)
                       Checkboxes{BadSensors(i)}.Visible='off'; 
                       h_axDec(BadSensors(i)).Visible='off';
                       h_axCln(BadSensors(i)).Visible='off';
                       h_axRaw(BadSensors(i)).Visible='off';
                       h_axRawDec(BadSensors(i)).Visible='off';
                    end
                    tempAx.YLim = [min(min(T(T>0)))-1,max(max(T(T>0))+1)];
                case 'Plot NaN (T = -999) data'
                    return
            end
        end 

    
                  
      %% Update loading box and close
        loading.Value =1;
        close(loading)

        
