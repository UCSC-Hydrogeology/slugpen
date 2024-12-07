%%% =======================================================================
%%  Purpose:
%       This function calculates tilt from acceleration data and corrects
%       depth if needed for SlugPen.
%%   Last edit:
%       07/19/2023 - Kristin Dickerson, UCSC
%%% =======================================================================

function [Tilt, G]...
            = TiltCalc(dataloaded)

%% Pull variables out of structures
    accz = dataloaded.ACCz; % Vectorized assignment
    z = dataloaded.Z; % Vectorized assignment

    %% Calculate Vertical Acceleration and Tilt
    G = accz * (2 / 256); % Vectorized calculation of G
    
    % Vertical Acceleration deviation from 1 G
    Gdev = abs(1 - abs(G)); % Vectorized deviation calculation
    
    % Tilt calculation: -1G = 0 deg; 0G = 90 deg
    Tilt = 90 * Gdev; % Linear mapping avoids intermediate steps
    
    %% Correct Depth (Vectorized)
    z(z < 0) = NaN; % Vectorized NaN assignment for invalid depths
    
    %% Save updated depth back to structure (if needed downstream)
    dataloaded.Z = z;