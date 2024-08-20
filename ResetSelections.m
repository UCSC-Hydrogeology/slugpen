%%% =======================================================================
%%  Purpose:
%     This function resets raw penetration in SlugPen so that  
%       no times are selected.
%%  Last edit
%     08/01/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function ResetSelections(H, dropdown_Select)

%% Delete lines on all plots
    delete(H.VerticalLines.Start_Eqm_Line);
    delete(H.VerticalLines.End_Eqm_Line);
    delete(H.VerticalLines.Start_Pen_Line);
    delete(H.VerticalLines.Start_Heat_Line);
    delete(H.VerticalLines.End_Pen_Line);
    delete(H.VerticalLines.SE_depth)
    delete(H.VerticalLines.EE_depth)
    delete(H.VerticalLines.SP_depth)
    delete(H.VerticalLines.HP_depth)
    delete(H.VerticalLines.EP_depth)
    delete(H.VerticalLines.SE_accel)
    delete(H.VerticalLines.EE_accel)
    delete(H.VerticalLines.SP_accel)
    delete(H.VerticalLines.HP_accel)
    delete(H.VerticalLines.EP_accel)
    delete(H.VerticalLines.SE_accel)
    delete(H.VerticalLines.EE_accel)
    delete(H.VerticalLines.SP_accel)
    delete(H.VerticalLines.HP_accel)
    delete(H.VerticalLines.EP_accel)
    delete(H.VerticalLines.EP_tilt)
    delete(H.VerticalLines.SE_tilt)
    delete(H.VerticalLines.EE_tilt)
    delete(H.VerticalLines.SP_tilt)
    delete(H.VerticalLines.HP_tilt)
    delete(H.VerticalLines.EP_tilt)

%% Reset drowdown and edit fields
    H.Selections.Value= dropdown_Select.Items(1);
    H.Selections.Start_Eqm.Value='';
    H.Selections.End_Eqm.Value='';
    H.Selections.Start_Pen.Value='';
    H.Selections.Start_Heat.Value='';
    H.Selections.End_Pen.Value='';