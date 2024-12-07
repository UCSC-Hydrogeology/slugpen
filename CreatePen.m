%%%% ======================================================================
%%  Purpose
%     This function creates the .pen file in SlugPen.
%
%%   Last edit
%      07/19/2023 by Kristin Dickerson, UCSC
%%% =======================================================================

function [HeatPulsePens,...
          DataType] ...
                = CreatePen(HeatPulsePens, ...
                    PulseData, ...
                    NoTherm, ...
                    figure_SetPenInfo, ...
                    figure_Main,...
                    DATA, ...
                    H, ...
                    PulsePower, ...
                    Cruisename, ...
                    StationName, ...
                    Penetration, ...
                    Latitude, ...
                    Longitude, ...
                    Datum,...
                    Depth_mean, ...
                    Tilt_mean, ...
                    LoggerID, ...
                    ProbeID, ...
                    NoActiveTherm)


%%  Initialize  

        % Define penetrations 
        % -------------------
        HeatPulsePens = [HeatPulsePens PulseData];
 
        % Define number of sensors
        % -------------------------
        NumThermTotal=NoTherm;

        % Ask user if they want to use raw, raw decimated, filtered, or
        % filtered decimate data to create the .pen file
        % ------------------------------------------------------------
        DataSelect = uiconfirm(figure_SetPenInfo, ['Would you like to ' ...
            'process the raw, filtered, decimated raw, or decimated ' ...
            'filtered data?'], 'Choose Data', 'Options',{'Raw', ...
            'Filtered', 'Raw decimated', 'Filtered decimated'});
 
        switch DataSelect
            case 'Raw' 
                 T      = DATA.Traw;
                 Tilt   = DATA.Tilt;
                 Depth  = DATA.Depth;
                 record = DATA.Record;
                 Time   = DATA.Time;
                 DataType = 1;
            case 'Filtered'
                 T      = DATA.Tcln;
                 Tilt   = DATA.Tilt;
                 Depth  = DATA.Depth;
                 record = DATA.Record;
                 Time   = DATA.Time;
                 DataType = 2;
            case 'Raw decimated'
                 T      = DATA.TdecRaw;
                 Tilt   = DATA.Tilt_dec;
                 Depth  = DATA.Depth_dec;
                 record = DATA.Record_dec;
                 Time   = DATA.Time_dec;
                 DataType = 3;
            case 'Filtered decimated'
                 T      = DATA.Tdec;
                 Tilt   = DATA.Tilt_dec;
                 Depth  = DATA.Depth_dec;
                 record = DATA.Record_dec;
                 Time   = DATA.Time_dec;
                 DataType = 4;
        end

         % Check year and correct if necessary 
         % -----------------------------------
            loggerStart = datevec(datetime(Time));
            loggerYear = loggerStart(1,1);
            if loggerYear < 2021
                loggerStart(:,1)=loggerStart(:,1)+2000;
                Time=datenum(loggerStart);
            end
          
%%   WRITE PENFILE  

            % Penfile Name
            % ------------
            cruise_name     = Cruisename;
            stn_name        = StationName;
            pen_name        = Penetration;
            penfile_name    = [cruise_name , '_', stn_name, '_', pen_name];

            PenetrationRecord_time = datetime(H.Selections.Start_Pen.Value); % START of PEN
            HeatPulseRecord_time   = datetime(H.Selections.Start_Heat.Value); % START of HEAT PULSE
            EndRecord_time         = datetime(H.Selections.End_Pen.Value); % END of entire RECORD

            BottomWaterStart_time   = datetime(H.Selections.Start_Eqm.Value); % START of bottom water CALIRBATION period
            BottomWaterEnd_time     = datetime(H.Selections.End_Eqm.Value); % END of bottom water CALIBRATION period

             % Penfile Records
             % ----------------
                % Variables for indexing, all decimated or regular time
                % intervals, depending which was selected by user
                
                % Records for calibrationa and penetration

                PenetrationRecords = record(Time>=PenetrationRecord_time & Time<=EndRecord_time); % records of penetration

                % Records for bottom water temps before penetration
                BWTime = find(Time<=PenetrationRecord_time); % Prior to penetration (only in bottom water)
                LastBWTime = BWTime(end); % final BW record before penetration
           
                % All records to be saved in PEN file
                PrePen = record(record<PenetrationRecords(1)); % Records prior to start of penetrtion
                PostPen = record(record>PenetrationRecords(end)); % Records following end of penetration
                PenfileRecords = [PrePen(end-5:end) PenetrationRecords PostPen(1:5)]; % All records to be used in PEN file, including 5 measurements before penetration and 5 measurements following end of penetration
                
             % Temperatures
             % ------------
                PenfileTemps       = T(1:NumThermTotal,record>=PenfileRecords(1) & record<=PenfileRecords(end)); % temps of each thermistor from start of penetration to end of record
                CalTemps           = T(1:NumThermTotal,Time>=BottomWaterStart_time & Time<=BottomWaterEnd_time); % temps of each thermistor from start of calibration to end of calibration
                BottomWaterRawData = T(1:NumThermTotal,LastBWTime); % temps of each thermistor right before penetration
                WaterSensorRawData = PenfileTemps(end, :); % temps of only bottom water thermistor ('T0')

                
                
            % Tapfile Records
            % ---------------
            TapfileRecords  = PenfileRecords;
            TapfileTilt     = abs(Tilt(record>=PenfileRecords(1) & record<=PenfileRecords(end)));
            TapfilePressure = Depth(record>=PenfileRecords(1) & record<=PenfileRecords(end));

            % Transpose all for writing
            % -------------------------
            BottomWaterRawData = BottomWaterRawData';
            PenfileRecords     = PenfileRecords';
            CalTemps           = CalTemps';
            PenfileTemps       = PenfileTemps';
            TapfileRecords     = TapfileRecords';
            TapfileTilt        = TapfileTilt';
            TapfilePressure    = TapfilePressure';
            WaterSensorRawData = WaterSensorRawData';

            % Offset times
            % -------------
                % Find closest record to penetration time
                foo = abs(PenetrationRecord_time-Time);
                a   = find(foo==min(foo));
                PenetrationRecord = record(a(1));
    
                % Find closest record to heat pulse time (if there was a pulse)
                if ~isnat(HeatPulseRecord_time)
                    foo = abs(HeatPulseRecord_time-Time);
                    a   = find(foo==min(foo));
                    HeatPulseRecord = record(a(1));
                else 
                    HeatPulseRecord=-999;
                end


            % Sequential records
            % -------------------
            PenfileRecords_sequential = 1:1:length(PenfileRecords);
            TapfileRecords_sequential = PenfileRecords_sequential;
            
            PenetrationRecord_sequential =  PenfileRecords_sequential(PenfileRecords==PenetrationRecord);
            HeatPulseRecord_sequential   = PenfileRecords_sequential(PenfileRecords==HeatPulseRecord);

            EndRecord_sequential         =  PenfileRecords_sequential(PenfileRecords==PenetrationRecords(end));

                % Set HeatPulseRecord_sequential to zero if no heat pulse
                if isempty(HeatPulseRecord_sequential)
                    HeatPulseRecord_sequential = -999;
                end
    
                % Get penetration number
                PenNum = str2double(extract(Penetration, digitsPattern));
    
                if ~isempty(PenNum)
                    PenPulsePower = num2str(PulsePower(PenNum));
                else
                    definepen=uiconfirm(figure_SetPenInfo, ['Please define ' ...
                        'penetration number.'], ...
                        'Must define metadata', 'Icon','warning', ...
                        'Interpreter','latex', 'Options','Okay');
                    waitfor(definepen)
                    return
                end

            

            % Write the .pen file
            % ---------------------
            PenfileName = [penfile_name,'.pen'];
            fido = fopen(['outputs/' PenfileName],'wt');
            if fido<1
                beep;
                uialert(figure_SetPenInfo, ['Unable to open ',PenfileName,' for writing. Check that you are in the correct directory.'], 'Error')
                return
            else

                % Header
                hdrstr1 = [StationName,' ',Penetration,' ','''', Cruisename,'''', ' ',Datum,];
                fprintf(fido,'%s\n',hdrstr1);
                fprintf(fido,'%6.1f %6.1f %6.0f %6.2f\n', str2double(Latitude), str2double(Longitude), str2double(Depth_mean), str2double(Tilt_mean));
                hdrstr3 = [LoggerID,'  ',ProbeID,'  ',num2str(NoActiveTherm), '  ', PenPulsePower];
                fprintf(fido,'%s\n',hdrstr3);
                fprintf(fido,'%6.0f\n',PenetrationRecord_sequential);
                fprintf(fido,'%6.0f \n',HeatPulseRecord_sequential);
    
                % Ouptut Format for Mean Calibration Temps
                Fmt = '%8.3f ' ;
                Fmt = repmat(Fmt,1,NumThermTotal);
                FmtBW = ['%6s',Fmt, '\n'];

                % Mean calibration temp data
                MeanCalTemps = mean(CalTemps, 1, 'omitnan');
                MeanCalTemps(isnan(MeanCalTemps))=-999;
                fprintf(fido,FmtBW, '  ', MeanCalTemps); 
                
                % Output Format for Thermistor Data
                Fmt = '%8.3f ' ;
                Fmt = repmat(Fmt,1,NumThermTotal);
                Fmt = ['%5.0f ',Fmt, '\n'];
    
                % Thermistor Data + tilt + depth (m)
                [nrows,~] = size(PenfileTemps);
                i=1;
                while i<=nrows
                    fprintf(fido,Fmt,PenfileRecords_sequential(i),PenfileTemps(i,:));
                    i=i+1;
                end
                fclose(fido);
            end

            % Write the .tap file
            % ---------------------
            TapfileName = [penfile_name,'.tap'];

            fido = fopen(['outputs/' TapfileName],'wt');
            if fido<1
                beep;
                disp(['Unable to open ',TapfileName,' for writing. Breaking.'])
                return
            else
    
                % Tilt and Pressure Data
                nrows = length(TapfileRecords);
                i=1;
                while i<=nrows
                    fprintf(fido,'%6.0f %6.2f %6.0f\n',TapfileRecords_sequential(i),TapfileTilt(i),TapfilePressure(i));
                    i=i+1;
                end
                fclose(fido);
            end

            % Save variables used in .pen file
            % --------------------------------------

            S_PENVAR = struct('BottomWaterRawData', BottomWaterRawData, 'AllRecords', ...
                PenfileRecords_sequential, 'CalibrationTemps', MeanCalTemps ,'AllSensorsRawData', PenfileTemps, 'WaterSensorRawData', WaterSensorRawData, ...
                'PenetrationRecord', PenetrationRecord_sequential, 'HeatPulseRecord', HeatPulseRecord_sequential, 'EndRecord', ...
                EndRecord_sequential);

            % Create data structure to store input and output values 
            S_PenHandles = struct('CruiseName',Cruisename,'StationName',StationName,...
                'Penetration',Penetration,'ProbeId',ProbeID,'Datum',Datum,...
                'LoggerId',LoggerID,'NoTherm',num2str(NoActiveTherm),'PulsePower', PenPulsePower, 'Latitude',Latitude,...
                'Longitude',Longitude,...
                'Tilt',Tilt_mean,'Depth',Depth_mean) ;

            
            % Save variables used in .tap file
            % --------------------------------------
            S_TAPVAR = struct('TAPRecord', TapfileRecords_sequential, 'Tilt', TapfileTilt, 'Depth', TapfilePressure);
            pause(1)
            
            % Save SlugPen.mat
            save(['outputs/' penfile_name], 'S_PENVAR', 'S_PenHandles', 'S_TAPVAR');
            MatfileName = [penfile_name '.mat'];

            % Save pen, tap, and mat files in the corresponding slugheat
            % input folder
            %cd outputs/
            %filenames = {PenfileName, TapfileName, MatfileName};
            %for i=1:length(filenames)
            %     copyfile(char(filenames(i)),'../../slugheat/inputs')
            %end
            %cd ..

            % Tell user if no calibration period was selected
            % ----------------------------------------------
            if all(MeanCalTemps==-999)
            conf = uiconfirm(figure_SetPenInfo, ['\bf No calibration period was selected ' newline newline ...
               'Temperature sensors will not be calibrated unless a ' newline ...
               'calibration period is chosen or mean calibration ' newline ...
               'temperatures for each sensor are manually input in ' ...
               '.pen file or .mat file by user'], ...
               'WARNING', 'Icon','warning', ...
               'Interpreter','latex', 'Options','Okay');
            waitfor(conf)
            end

            figure_SetPenInfo.Visible='off';

            % Tell user where the outputs are saved
            % ----------------------------------------
            uialert(figure_Main, ['Output files created for penetration: ' newline newline...
               penfile_name '.pen' newline ...
               penfile_name '.tap' newline ...
               penfile_name '.mat' newline newline...
               'Files saved in in your current path "outputs" subfolder ' pwd '/outputs/'], ...
               'Success', 'Icon','success')




 