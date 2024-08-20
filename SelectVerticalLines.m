
%%% =======================================================================
%%  Purpose: 
%       This function allows the user to select timing of events in SlugPen.
%% Last edit:
%       01/22/2024 by Kristin Dickerson, UCSC
%%% =======================================================================

function [lines] = SelectVerticalLines(H, ...
            figure_Main, ...
            dropdown_Select)

%% Initialize
        zoom(H.Axes.Raw, 'off')
        pan(H.Axes.Raw, 'off')
        enableDefaultInteractivity(H.Axes.Raw);
        x = dropdown_Select.Value;

        % Ensure user is in compatible MATLAB version
        % --------------------------------------------
        if verLessThan('Matlab', '9.9')
            uialert(figure_Main, ['This function is not supported in Matlab ' ...
                'versions prior to Matlab r2020b.'], 'Error Matlab Version')
            return
        end

%% Allow user selection of time using crosshair selector (ginput function)        

         % Set up figure handle visibility
         % --------------------------------
        fhv = figure_Main.HandleVisibility;        % Current status
        figure_Main.HandleVisibility = 'callback'; % Temp change (or, 'on') 
        set(0, 'CurrentFigure', figure_Main)       % Make fig current
        
        % Run ginput to allow selection
        % -----------------------------
        [selTime,~] = ginput(1);                     
    
        % Return figure to original state
        % --------------------------------
        figure_Main.HandleVisibility = fhv;       
      
%% Plot selection lines on each plot (temp, depth, and tilt axes)   
        
        % Plot
        % ----
        depthLine = xline(H.Axes.Depth, selTime, '-.k', 'LineWidth', 2);
        tiltLine = xline(H.Axes.Tilt, selTime, '-.k', 'LineWidth', 2);
        accelLine = xline(H.Axes.Accel, selTime, '-.k', 'LineWidth', 2);
        VerticalLine = xline(H.Axes.Raw, selTime, '-.k', x, 'LineWidth', 2, ...
            "LabelHorizontalAlignment","center","FontSize",24, "FontWeight","bold");
        
        % Set axes
        % --------
        VerticalLine.Label = x;
        VerticalLine.LabelHorizontalAlignment = 'center';
        drawnow;

        % Convert to a datetime
        % -----------------------
        selTime = num2ruler(selTime,H.Axes.Raw.XAxis);

        % Save lines
        % -----------
        SE_depthLine = H.VerticalLines.SE_depth;
        SE_accelLine = H.VerticalLines.SE_accel;
        SE_tiltLine  = H.VerticalLines.SE_tilt;
        EE_depthLine = H.VerticalLines.EE_depth;
        EE_accelLine = H.VerticalLines.EE_accel;
        EE_tiltLine  = H.VerticalLines.EE_tilt;
        SP_depthLine = H.VerticalLines.SP_depth;
        SP_accelLine = H.VerticalLines.SP_accel;
        SP_tiltLine  = H.VerticalLines.SP_tilt;
        HP_depthLine = H.VerticalLines.HP_depth;
        HP_accelLine = H.VerticalLines.HP_accel;
        HP_tiltLine  = H.VerticalLines.HP_tilt;
        EP_depthLine = H.VerticalLines.EP_depth;
        EP_accelLine = H.VerticalLines.EP_accel;
        EP_tiltLine  = H.VerticalLines.EP_tilt;
       
%% Label lines based on which event is selected in the drop down menu
        switch x
            case 'Start Calibration Period (C1)'
                H.Selections.Start_Eqm.Value = string(selTime);
                delete(H.VerticalLines.Start_Eqm_Line)
                H.VerticalLines.Start_Eqm_Line = VerticalLine;
                delete(SE_depthLine)
                SE_depthLine = depthLine;
                delete(SE_accelLine)
                SE_accelLine = accelLine;
                delete(SE_tiltLine)
                SE_tiltLine = tiltLine;
            case 'End Calibration Period (C2)'
                H.Selections.End_Eqm.Value = string(selTime);
                delete(H.VerticalLines.End_Eqm_Line)
                H.VerticalLines.End_Eqm_Line = VerticalLine;
                delete(EE_depthLine)
                EE_depthLine = depthLine;
                delete(EE_accelLine)
                EE_accelLine = accelLine;
                delete(EE_tiltLine)
                EE_tiltLine = tiltLine;
            case 'Start Penetration (P)'
                H.Selections.Start_Pen.Value = string(selTime);
                delete(H.VerticalLines.Start_Pen_Line)
                H.VerticalLines.Start_Pen_Line = VerticalLine;
                delete(SP_depthLine)
                SP_depthLine = depthLine;
                delete(SP_accelLine)
                SP_accelLine = accelLine;
                delete(SP_tiltLine)
                SP_tiltLine = tiltLine;
            case 'Heat Pulse (H)'
                H.Selections.Start_Heat.Value = string(selTime);
                delete(H.VerticalLines.Start_Heat_Line)
                H.VerticalLines.Start_Heat_Line = VerticalLine;
                delete(HP_depthLine)
                HP_depthLine = depthLine;
                delete(HP_accelLine)
                HP_accelLine = accelLine;
                delete(HP_tiltLine)
                HP_tiltLine = tiltLine;
            case 'End Penetration (E)'
                H.Selections.End_Pen.Value = string(selTime);
                delete(H.VerticalLines.End_Pen_Line)
                H.VerticalLines.End_Pen_Line= VerticalLine;
                delete(EP_depthLine)
                EP_depthLine = depthLine;
                delete(EP_accelLine)
                EP_accelLine = accelLine;
                delete(EP_tiltLine)
                EP_tiltLine = tiltLine;
        end
        zoom(H.Axes.Raw, 'on')
        pan(H.Axes.Raw, 'on')
    
%% Structure for lines
    lines = struct('Start_Eqm_Line',H.VerticalLines.Start_Eqm_Line,...
        'End_Eqm_Line', H.VerticalLines.End_Eqm_Line, 'Start_Pen_Line', ...
        H.VerticalLines.Start_Pen_Line, 'Start_Heat_Line', ...
        H.VerticalLines.Start_Heat_Line, 'End_Pen_Line', H.VerticalLines.End_Pen_Line, ...
        'SE_depth',SE_depthLine, 'SE_accel',SE_accelLine, 'SE_tilt',SE_tiltLine, ...
        'EE_depth',EE_depthLine, 'EE_accel',EE_accelLine, ...
        'EE_tilt',EE_tiltLine,'SP_depth',SP_depthLine,'SP_accel',SP_accelLine, ...
        'SP_tilt',SP_tiltLine,'HP_depth',HP_depthLine, ...
        'HP_accel',HP_accelLine,'HP_tilt',HP_tiltLine,'EP_depth',EP_depthLine, ...
        'EP_accel',EP_accelLine,'EP_tilt',EP_tiltLine);





