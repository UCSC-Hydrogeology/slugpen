%%% ==============================================================================
%   Purpose: 
%     This function turns sensor lines on and off based on sensor legend
%     checkbox values. When turning the line back ON, only the decimated
%     line will be replotted. This needs to be revisited in future if you
%     want to be able to replot raw, clean, and decimated plots.
%%% ==============================================================================

function cbSensValueChange(~, src, Dec, Clean, Raw, Rawdec, CleanCb, DecCb, RawCb, RawDecCb)
            if src.Value == 0               % if checkbox is unchecked
                Dec.Visible = 'off';
                Clean.Visible = 'off';
                Raw.Visible = 'off';
                Rawdec.Visible = 'off';
            else

                Raw.Visible = 'on';


            end
            
