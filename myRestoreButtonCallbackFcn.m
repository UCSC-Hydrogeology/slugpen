% Responds to pressing the restore button in the axes' toolbar. 
% originalCallback is a function handle to the original callback 
%   function for this button. 

%%%% ======================================================================
%%  Purpose
%     Responds to pressing the restore button in the axes' toolbar. 
%     originalCallback is a function handle to the original callback 
%     function for this button. 
%
%%   Last edit
%      07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function myRestoreButtonCallbackFcn(hobj, event, originalCallback,H)

% Evaluate original callback
    originalCallback(hobj,event) 

% Update axes  
    H.Axes.Raw.XLim = H.Axis_Limits.AxLims_Raw;