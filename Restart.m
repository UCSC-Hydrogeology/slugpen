%%% =======================================================================
%%  Purpose:
%       This function resets SlugPen to its default.
%%  Last edit
%       07/19/2023 by Kristin Dickerson, UCS
%%% =======================================================================

function [All_Zeros, ...
          All_NaN,...
          Traw,...
          Tcln,...
          Pitch,...
          Roll,...
          G,...
          Tilt,...
          Time,...
          AccelPlot,...
          TiltPlot,...
          DepthPlot,...
          penetrationInfo,...
          HeatPulsePens,...
          PulseData,...
          selTime_Vector,...
          DataType,...
          Params]...
            = Restart(figure_Main,...
            H, ...
            grid_PlotWindow, ...
            checkbox_UseHP,...
            PulseData, ...
            label_currentpathfull, ...
            dropdown_Select)


%% Constants and Initializations
    
    if isempty(H.Fileinfo.No_Thermistors)
        NoTherm=[];
    else
        NoTherm = H.Fileinfo.No_Thermistors.Value; % Get number of thermistors
    end
    All_Zeros = zeros(1, NoTherm);
    All_NaN = nan(1, NoTherm); % Use NaN directly instead of division

%% RESET CONTROLS
    H.Error=0;
     H.Exe_Controls.Delete.Enable='off';    
    controlsToDisable = {'Delete', 'RawPlot', 'CleanPlot', 'DecimatedPlot', 'RawDecimatedPlot'};
    for control = controlsToDisable
        H.Plot_Controls.(control{1}).Enable = 'off';
    end
     
 %% CLEAR AXES 
    axesFields = {'Raw', 'Depth', 'Tilt', 'Accel'};
    for field = axesFields
        cla(H.Axes.(field{1}), 'reset');
        H.Axes.(field{1}).XLim = datenum([datetime('yesterday') datetime('today')]);
    end

     % Single-time dummy plot for all axes
    x = datetime('yesterday'):minutes(60):datetime('today');
    y = linspace(0.04, 1, numel(x));
    for field = axesFields
        plot(H.Axes.(field{1}), x, y, 'Visible', 'off');
    end

     AccelPlot = 0;
     TiltPlot = 1;
     DepthPlot = 1;
     
%% RESET CHECKBOXES
    checkboxFields = {'Accel_check', 'Depth_check', 'Tilt_check', ...
                      'Raw_check', 'Rawdec_check', 'Filtered_check', ...
                      'Filtereddec_check'};
    checkboxValues = [0, 1, 1, 1, 1, 0, 0];
    for i = 1:numel(checkboxFields)
        H.PlotCheckboxes.(checkboxFields{i}).Value = checkboxValues(i);
    end
     
    %% RESET MAIN GRID
    grid_PlotWindow.RowHeight = {'0.2x', '0.27x', '5x', '1x', '0x', '1x', '0.2x', '0.2x'};

     % Vertical line begins at x = 0
     enableDefaultInteractivity(H.Axes.Raw);
    
%% RESET COMPONENTS
    selectionsFields = {'Raw', 'Depth', 'Tilt', 'Accel', ...
                        'Raw_Text', 'Depth_Text', 'Tilt_Text', 'Accel_Text'};
    for field = selectionsFields
        H.Selections_Lines.(field{1}) = All_Zeros;
    end
     
    % Reset Export/Write Buttons
    H.Exe_Controls.Export_Penfile.Enable = 'off';

     %% RESET FILE INFORMATION
    fileInfoFields = {'Filename', 'Start_Date', 'Start_Time', 'End_Time', ...
                      'Log_Interval', 'Pulse_Length', 'Pulse_Power', 'Decay_Time'};
    for field = fileInfoFields
        H.Fileinfo.field{1}.Value = '';
    end

     %% RESET SELECTION INFO
    H.Selections.Select.Enable = 'on';
    dropdownFirstItem = dropdown_Select.Items{1};
    H.Selections.Value = dropdownFirstItem;
    selectionFields = {'Start_Eqm', 'End_Eqm', 'Start_Pen', 'Start_Heat', 'End_Pen'};
    for field = selectionFields
        H.Selections.(field{1}).Value = '';
    end

  %% RESET PENETRATION INFORMATION

    penetrationInfo = strings(1, 11);

     % HEAT PULSE DATA
    HeatPulsePens = [];
    PulseData = checkbox_UseHP.Value;
  
%% DEFAULT SCREEN
    figure_Main.WindowState = 'maximized';
      

%% INITIALIZE VECTORS
      
      %% INITIALIZE VECTORS
    H.S_Selections_Lines = struct('Raw', All_NaN, 'Depth', All_NaN, ...
                                  'Tilt', All_NaN, 'Raw_Text', All_NaN, ...
                                  'Depth_Text', All_NaN, 'Tilt_Text', All_NaN);
    H.S_Plot_Controls.Legend = All_NaN;
    H.S_AxLims = struct('AxLims_Raw', All_NaN, 'AxLims_Depth', All_NaN, ...
                        'AxLims_Tilt', All_NaN, 'AxLims_Acc', All_NaN);
    selTime_Vector = All_NaN;

    % Initialize scalar variables
    [Traw, Tcln, Pitch, Roll, G, Tilt, Time] = deal(0);
    [H.Error, H.Crop] = deal(0);
    DataType = 4;
    Params = {'text'};

    %% SET CURRENT PATH
    CurrentPath = pwd;
    label_currentpathfull.Text = ['...' CurrentPath(end-20:end)];

        H.Axes.Accel.Visible = 'off';
end

         