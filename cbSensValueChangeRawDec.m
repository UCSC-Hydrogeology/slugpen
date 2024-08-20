%%% =======================================================================
%% Purpose: 
%     This function turns sensor lines on and off based on sensor legend
%     checkbox values. When turning the line back ON, only the raw, 
%     decimated line will be replotted. 
%% Last edit:
%     03/12/2024 by Kristin Dickerson, UCSC
%%% =======================================================================

function cbSensValueChangeRawDec(~, src, Dec, Clean, Raw, Rawdec, ~, ~, ~, ~)
% if checkbox is unchecked, turn off plots                
            if src.Value == 0               
                Dec.Visible = 'off';
                Clean.Visible = 'off';
                Raw.Visible = 'off';
                Rawdec.Visible = 'off';
% if checkbox is checked, turn on only raw, decimated plot                    
            else

                Rawdec.Visible = 'on';

            end
            