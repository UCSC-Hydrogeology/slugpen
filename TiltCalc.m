%%% =======================================================================
%%  Purpose:
%       This function calculates tilt from acceleration data and corrects
%       depth if needed for SlugPen.
%%   Last edit:
%       07/19/2023 - Kristin Dickerson, UCSC
%%% =======================================================================

function [Tilt, G]...
            = TiltCalc(dataloaded, ...
                    parsedtiming)

%% Pull variables out of structures
         accz   = dataloaded.ACCz;
         z      = dataloaded.Z;

         yr    = parsedtiming.YR;

 %% CALCULATE TILT 
            
            % Vertical Acceleration: -256 --> +256 == -2g --> + 2g
            G = accz*(2/256);

            % Vertical Acceleration deviation from 1 G
            Gdev = abs(1-abs(G));

            % Calculate Tilt : -1G = 0 deg; 0G = 90 deg
            vert = 0; % when probe is vertical
            hor = 90; % when probe is horizontal
            devMin = 0; % minumum deviation from vertical
            devMax = 1; % maximum deviation from vertical
            
            Tilt = vert-((Gdev-devMin))*(vert-hor)./(devMax-devMin);

%% CORRECT DEPTH, if needed
        a=find(z<0); % changed Z's that are less than 0 to NaN
        if ~isempty(a)
            z(a)=NaN;
        end    
         
%% Save variables back in structures
         dataloaded.ACCz = accz;
         dataloaded.Z    = z;

         parsedtiming.YR = yr;