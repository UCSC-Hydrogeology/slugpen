%%% ==============================================================================
%%  Purpose: 
%     This function allows user to select timing of events using ROI
%%  Last edit:
%     07/19/2023 by Kristin Dickerson, UCSC
%%% ==============================================================================

function [rawData,...
          clnData,...
          decData,...
          decRaw]...
            = SelectROI(H,DATA)
    
     NoTherm = H.Fileinfo.No_Thermistors.Value;
     
     % Turn off zoom and pan
     pan(H.Axes.Raw, "off");
     zoom(H.Axes.Raw, "off");
     
     % Extract Time
     Time     = DATA.Time;       % For Raw and Cln
     Time_dec = DATA.Time_dec;   % For Decimated Raw and Cln
     
     % Extract temp data
     Traw = DATA.Traw;
     Tcln = DATA.Tcln;
     Tdec = DATA.Tdec;
     TdecRaw = DATA.TdecRaw;
     
     % Combine time and temp data
     rawData = [Time;Traw];
     clnData = [Time; Tcln];
     decData = [Time_dec; Tdec];
     decRaw = [Time_dec; TdecRaw];
     disableDefaultInteractivity(axis_raw)
     roiSelect = drawfreehand(axis_raw);
     roiSelect.InteractionsAllowed = 'none';
     drawnow
     
     % Selected xdata and ydata %
     i=1;
     rawSelected=zeros(5174, NoTherm);
     clnSelected=zeros(5174, NoTherm);
     decSelected=zeros(518, NoTherm);

     while i<=NoTherm
         rawSelected(:,i) = inROI(roiSelect, Time, Traw(i,:));
         i=i+1;
     end
     
     rawSelected=rawSelected(i,:);