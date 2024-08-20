%%% =======================================================================
%%  Purpose:
%     This function aligns the inner position left edge and width (x axis)
%     of all 4 plots of SlugPen
%%  Last edit
%     07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function xAlign(axis_raw,...
                axis_depth, ...
                axis_accel, ...
                axis_tilt)

% Pause so code can catch up
            pause(.2);
  
%% Align all axes      
            axis_accel.InnerPosition(1) = axis_depth.InnerPosition(1);
            axis_accel.InnerPosition(3) = axis_depth.InnerPosition(3);
            axis_raw.InnerPosition(1) = axis_depth.InnerPosition(1);
            axis_raw.InnerPosition(3) = axis_depth.InnerPosition(3);
            axis_tilt.InnerPosition(1) = axis_depth.InnerPosition(1);
            axis_tilt.InnerPosition(3) = axis_depth.InnerPosition(3);

% Pause so code can catch up
            pause(.2);
            drawnow;