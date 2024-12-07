%%% =======================================================================
%%  Purpose
%       This function filters the raw temperature data. This includes 
%       filtering outliers, smoothing, and decimating data. Filtering 
%       method by Mike Hutnak, ROQ.
%%   Last edit:
%       07/19/2023 - Kristin Dickerson, UCSC
%%% =======================================================================

function [Tcln, ...
          decimateddata,...
          Record] ...
            = Filter(figure_Main, ...
            edit_MedianFilter,...
            edit_MeanFilter,...
            edit_DecimateFilter,...
            datauniquetimes, ...
            Traw, ...
            NoTherm)


%% Get filtering windows from user input in GUI

            wlmedian = edit_MedianFilter.Value;
            wlmean = edit_MeanFilter.Value;
            wldec = edit_DecimateFilter.Value;

             if rem(wlmean,2)~=0
                uialert(figure_Main, 'Mean window length must be even number',...
                'Filter Error');
                return
             end

%% Pull data from structure
            timeNumU   =  datauniquetimes.TIMENUMU ;
            pwrU       =  datauniquetimes.PWRU ;
            tiltU      =  datauniquetimes.TILTU ;
            gU         =  datauniquetimes.GU ;
            zU         =  datauniquetimes.ZU ;
            TrawU      =  datauniquetimes.TRAWU ;
       
            
%% RAW DATA CLEANING
         Tallf1 = Traw;
         Tallf1(Tallf1 <= 0) = NaN; % Clean invalid data points (<= 0)

%% FILTER DATA 
            % Despike and smooth the data for each sensor in one step
            Tallf2 = Tallf1; % Start with the raw data

            % Apply median filter to each sensor across columns
            Tallf2 = movmedian(Tallf2, wlmedian, 2); 
            % Apply mean filter to the data after median filtering
            Tallf2 = rbmmed(Tallf2, wlmean); 
            
            % Assign filtered data to Tcln
            Tcln = Tallf2;

%% DECIMATE DATA
            n = 1:wldec:length(Traw(1,:)); % Indices for decimation
            TdecRaw = Traw(:, n); % Decimated raw data
            Tdec = Tcln(:, n); % Decimated filtered data
            timeDec = timeNumU(n); % Decimated time
            
%% Finalize
            % Create Record Counter
            Record = 1:length(timeNumU);
            
            % Prepare decimated data
            Record_dec = Record(n);
            G_dec = gU(n);
            Tilt_dec = tiltU(n);
            Power_dec = pwrU(n);
            Depth_dec = zU(n);
            
            % Save decimated data in structure
            decimateddata = struct('TDecRaw', TdecRaw, 'TDecFiltered', Tdec, ...
                'TimeDec', timeDec, 'RecordDec', Record_dec, 'GDec', G_dec, ...
                'TiltDec', Tilt_dec, 'PowerDec', Power_dec, 'DepthDec', Depth_dec);
            
end

