%%% =======================================================================
%%  Purpose:
%       This function creates datatimes for data loaded into SlugPen as well 
%       as takes only unique times for further processing
%%   Last edit:
%        07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function [datauniquetimes] ...
            = GetDateTime(dataloaded,parsedtiming, ...
                        Tilt,G) 

%% Pull variables out of structures
         accx   = dataloaded.ACCx;
         accy   = dataloaded.ACCy;
         accz   = dataloaded.ACCz;
         pwr    = dataloaded.PWR;
         z      = dataloaded.Z;
         Traw   = dataloaded.TRAW;   
         yr    = parsedtiming.YR;
         mo    = parsedtiming.MO;
         dy    = parsedtiming.DY;
         hr    = parsedtiming.HR;
         mn    = parsedtiming.MN;
         sc    = parsedtiming.SC;

 %% Correct year, if needed
        if yr(1)<2000
            yr=yr+2000;
        end

        timeNum  = datetime(yr,mo,dy,hr,mn,sc); 
        [~,indU] = unique(timeNum);  
        timeNumU = timeNum(indU);
        accxU    = accx(indU);
        accyU    = accy(indU);
        acczU    = accz(indU);
        pwrU     = pwr(indU);
        tiltU    = Tilt(indU);
        gU       = G(indU);
        zU       = z(indU);
        TrawU    = Traw(:,indU);


%% Save unique data in a structure for access
        datauniquetimes = struct('TIMENUMU', timeNumU, 'ACCxU', accxU, ...
            'ACCyU', accyU, 'ACCzU', acczU, 'PWRU', pwrU, 'TILTU',...
            tiltU, 'GU',gU, 'ZU',zU, 'TRAWU',TrawU);

