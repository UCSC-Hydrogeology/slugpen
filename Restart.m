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

 %% RESET CONTROLS
     H.Error=0;
     H.Exe_Controls.Delete.Enable='off';
     H.Plot_Controls.RawPlot.Enable='off';
     H.Plot_Controls.CleanPlot.Enable='off';
     H.Plot_Controls.DecimatedPlot.Enable='off';
     H.Plot_Controls.RawDecimatedPlot.Enable='off';
     
     % Initialize Selection Vectors for 22 if No_Thermistors hasn't been
     % populated yet

      % Initial zeros or NaN
      NoTherm = 11;
      All_Zeros = zeros(1,NoTherm);
      All_NaN = All_Zeros/0;
     
 %% CLEAR AXES 
     cla(H.Axes.Raw, 'reset');
     cla(H.Axes.Depth, 'reset');
     cla(H.Axes.Tilt, 'reset');
     
     % Make all axes visible
     H.Axes.Accel.Visible='off';
     H.Axes.Tilt.Visible='on';
     H.Axes.Depth.Visible='on';
     H.Axes.Raw.Visible='on';   

     x=datetime('yesterday'):minutes(60): datetime('today');
     y=1/25:1/25: 1;

     plot(H.Axes.Raw,x,y, 'Visible','off')
     plot(H.Axes.Accel,x,y,'Visible','off')
     plot(H.Axes.Depth,x,y,'Visible','off')
     plot(H.Axes.Tilt,x,y,'Visible','off')

     H.Axes.Raw.XLim = [datetime('yesterday') datetime('today')]; 
     H.Axes.Accel.XLim = [datetime('yesterday') datetime('today')];
     H.Axes.Tilt.XLim = [datetime('yesterday') datetime('today')];
     H.Axes.Depth.XLim = [datetime('yesterday') datetime('today')];

     AccelPlot = 0;
     TiltPlot = 1;
     DepthPlot = 1;
     
     % Turn on plot control checkboxes
     H.PlotCheckboxes.Accel_check.Value=0;
     H.PlotCheckboxes.Depth_check.Value=1;
     H.PlotCheckboxes.Tilt_check.Value=1;
     H.PlotCheckboxes.Raw_check.Value=1;
     H.PlotCheckboxes.Rawdec_check.Value=1;
     H.PlotCheckboxes.Filtered_check.Value=0;
     H.PlotCheckboxes.Filtereddec_check.Value=0;   
     
     % Reset main grid
     grid_PlotWindow.RowHeight={'0.2x','0.27x','5x','1x','0x','1x','0.2x','0.2x'};   

     % Vertical line begins at x = 0
     enableDefaultInteractivity(H.Axes.Raw);
    
  %% RESET ALL COMPONENTS 
         H.Selections_Lines.Raw                   = All_Zeros;
         H.Selections_Lines.Depth                 = All_Zeros;
         H.Selections_Lines.Tilt                  = All_Zeros;
         H.Selections_Lines.Accel                 = All_Zeros;
         H.Selections_Lines.Raw_Text              = All_Zeros;
         H.Selections_Lines.Depth_Text            = All_Zeros;
         H.Selections_Lines.Tilt_Text             = All_Zeros;
         H.Selections_Lines.Accel_Text            = All_Zeros;   
     
     % Reset Export/Write Buttons
         H.Exe_Controls.Export_Penfile.Enable='off';

     % Reset Datafile Information
         H.Fileinfo.Filename.Value = '';
         H.Fileinfo.Start_Date.Value='';
         H.Fileinfo.Start_Time.Value='';
         H.Fileinfo.End_Time.Value='';
         H.Fileinfo.Log_Interval.Value='';
         H.Fileinfo.Pulse_Length.Value='';
         H.Fileinfo.Pulse_Power.Value='';
         H.Fileinfo.Decay_Time.Value='';   

     % Reset Selection info
         H.Selections.Select.Enable='on';
         H.Selections.Value=dropdown_Select.Items(1);
         H.Selections.Start_Eqm.Value='';
         H.Selections.End_Eqm.Value='';
         H.Selections.Start_Pen.Value='';
         H.Selections.Start_Heat.Value='';
         H.Selections.End_Pen.Value='';

  %% RESET PENETRATION INFORMATION

     % Penetration Info %
     penetrationInfo = ["", "", "", "", "", "", "", "" ...
         "", "", ""];

      pause(2)

      % Initialize array with heat pulse data
      HeatPulsePens = [];
      Pulse = checkbox_UseHP.Value;
      if Pulse == 1
          PulseData = 1;
      elseif Pulse == 0
          PulseData = 0;
      end

  
%% DEFAULT SCREEN
      % Make main figure full screen %
      figure_Main.WindowState = 'maximized';
      

%% INITIALIZE VECTORS
      
      % Set Initial NaN conditions %
      H.S_Selections_Lines.Raw = All_NaN;
      H.S_Selections_Lines.Depth = All_NaN;
      H.S_Selections_Lines.Tilt = All_NaN;
      H.S_Selections_Lines.Raw_Text = All_NaN;
      H.S_Selections_Lines.Depth_Text = All_NaN;
      H.S_Selections_Lines.Tilt_Text = All_NaN;

      H.S_Plot_Controls.Legend = All_NaN;

      H.S_AxLims.AxLims_Raw   = All_NaN;
      H.S_AxLims.AxLims_Depth = All_NaN;
      H.S_AxLims.AxLims_Tilt  = All_NaN;
      H.S_AxLims.AxLims_Acc   = All_NaN;

      selTime_Vector = All_NaN;

      % Set Initial Zero Conditions %
      Traw = 0;
      Tcln = 0;
      Pitch = 0;
      Roll = 0;
      G = 0;
      Tilt = 0;
      Time = 0;
      H.Error = 0;
      H.Crop = 0;
      DataType = 4;

      % Initial Parameter labels %
       Params = {'text'}; 

       % Set Current path
       % -----------------
       CurrentPath = pwd;
       label_currentpathfull.Text = ['...' CurrentPath(end-20:end)];
        

         