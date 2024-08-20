%%% =================================================================
%%  Purpose:
%       This function resets raw penetration in SlugPen so that  
%       penetration times are not selected, but calibration times remain.
%%  Last edit
%       01/22/2024 by Kristin Dickerson, UCSC
%%% =================================================================

function ResetPen(H)

%% Delete lines on all plots
         delete(H.VerticalLines.Start_Pen_Line);
         delete(H.VerticalLines.Start_Heat_Line);
         delete(H.VerticalLines.End_Pen_Line);

         delete(H.VerticalLines.SP_depth)
         delete(H.VerticalLines.HP_depth)
         delete(H.VerticalLines.EP_depth)

         delete(H.VerticalLines.SP_accel)
         delete(H.VerticalLines.HP_accel)
         delete(H.VerticalLines.EP_accel)

         delete(H.VerticalLines.SP_accel)
         delete(H.VerticalLines.HP_accel)
         delete(H.VerticalLines.EP_accel)
         delete(H.VerticalLines.EP_tilt)

         delete(H.VerticalLines.SP_tilt)
         delete(H.VerticalLines.HP_tilt)
         delete(H.VerticalLines.EP_tilt)

%% Reset drowdown and edit fields
         H.Selections.Value=dropdown_Select.Items(1);

         H.Selections.Start_Pen.Value='';
         H.Selections.Start_Heat.Value='';
         H.Selections.End_Pen.Value='';