%%% =======================================================================
%%  Purpose:
%       This function loads data on Unix systems for SlugPen.
%%  Last edit:
%       07/19/2023 - Kristin Dickerson, UCSC
%%% =======================================================================

function [Parameters, ...
          dataloaded,...
          parsedtiming,...
          HSTART,...
          RPT] ...
                = LoadUnix(figure_Main, ...
                    NoTherm, ...
                    f)

%% Initialize
     
     % Get file name
     [~,fn,~]=fileparts(f);
        
%% Scan using UNIX awk
                    
           % Delete temp file
           if exist('RawData.tmp','file')
               delete('RawData.tmp');
           end

           % Delete temp file
           if exist('RawDataTemp.tmp','file')
               delete('RawDataTemp.tmp');
           end

            % Delete temp file
           if exist('HFParams.tmp','file')
               delete('HFParams.tmp');
           end

           % Delete temp file
           if exist([fn,'_HReportTemp.dat'],'file')
               delete([fn,'_HReportTemp.dat']);
           end
           
           
           % Grab all header lines
           eval(['!awk ''/@/'' ',f,' > ',fn,'_Header.dat'])
           % Grab all lines with Data
           eval(['!awk ''!/^@/'' ',f,' > RawDataTemp.tmp'])
           rtemp = fullfile(pwd,'RawDataTemp.tmp');
           fid = fopen(rtemp,'r');
           eval(['!awk ''length > 180'' ',rtemp,' > RawData.tmp'])

           % Remove reports with negative power
           % Grab all Acomms Report lines
           eval(['!awk ''/RPT/'' ',f,' > ',fn,'_HReportTemp.dat'])
           % Define the input and output filenames
                inputFileName = [fn,'_HReportTemp.dat']; % Replace with your actual input file name
                outputFileName = [fn,'_HReport.dat']; % Replace with your desired output file name
                
                % Open the input file for reading
                fileID = fopen(inputFileName, 'r');
                if fileID == -1
                    error('Cannot open the input file.');
                end
                
                % Read all lines from the file into a cell array
                lines = {};
                while ~feof(fileID)
                    line = fgetl(fileID);
                    if ischar(line)
                        lines{end+1} = line;
                    end
                end
                fclose(fileID);
                
                % Filter lines that do not contain 'J-'
                filteredLines = lines(~contains(lines, 'J-'));
                
                % Open the output file for writing
                fileID = fopen(outputFileName, 'w');
                if fileID == -1
                    error('Cannot open the output file.');
                end
                
                % Write filtered lines to the output file
                for i = 1:numel(filteredLines)
                    fprintf(fileID, '%s\n', filteredLines{i});
                end
                fclose(fileID);

           % Grab all Heat Pulse Report lines
           eval(['!awk ''/HSTART/'' ',f,' > ',fn,'_HPulse.dat'])

           % Grab all Changes Report lines
           eval(['!awk ''/STATECHANGE/'' ',f,' > ',fn,'_ChangesReport.dat'])
           
           % Prep textscan for variable number of thermistors
           n  = repmat('%n ',1,NoTherm); %Thermistors
           s       = repmat('%s ',1,11);     % Time, Acc, depth
           fid = fopen('RawData.tmp','r');
           C   = textscan(fid,[s n]);
           fclose(fid);
      
           % Initialize Matrix AllTdata (All Temperature Data)
           l    = length(C{1});
           Traw = NaN*zeros(NoTherm+1,l);
      
           % Assign parameters
           ymd    = str2num(char(strrep(C{3},'-',' '))); 
           hms    = str2num(char(strrep(C{4},':',' ')));
           accx   = str2num(char(strrep(C{6},'X','')));
           accy   = str2num(char(strrep(C{7},'Y','')));
           accz   = str2num(char(strrep(C{8},'Z','')));
           pwr    = str2num(char(strrep(C{9},'V','')));
           z      = str2num(char(strrep(C{10},'D','')));
      
           for i=1:NoTherm
               foo = C{11+i}';
               Traw(i,1:length(foo))=foo;
           end
      
           % Restructure Traw with T1 in row 1, incrementing to Twater
           Traw(NoTherm+1,:)=Traw(1,:);
           Traw=Traw(2:end,:);

           % Save data in structure for access
            dataloaded = struct('YMD',ymd,'HMS', hms,'ACCx', accx,...
           'ACCy',accy,'ACCz',accz,'PWR',pwr,'Z',z, 'TRAW',Traw);
     
       %% Error if number of temperature measurements ~= time
       if length(Traw)~=length(hms)
           errorstr='Timing Error';
               dlgname = 'Timing Error';
           uiconfirm(figure_Main, errorstr, dlgname, ...
               'Options', {'OK','Cancel'});
           return
       end  
       
       %% Parse Timing
       yr    = ymd(:,1);
       mo    = ymd(:,2);
       dy    = ymd(:,3);
       hr    = hms(:,1);
       mn    = hms(:,2);
       sc    = hms(:,3);

       % Save data in structure for access
       parsedtiming = struct('YR',yr,'MO', mo,'DY', dy,...
               'HR',hr,'MN',mn,'SC',sc);

       %% Read in heat pulse and report data
       HSTART = {};
       fid = fopen([fn,'_HPulse.dat'],'r');
       i=1;
       while i>0
           line = fgetl(fid);
           if line == -1
               fclose(fid);
               break
           else
               HSTART{i}=strtrim(line); %#ok<AGROW>
               i=i+1;
           end
       end
       RPT = {};
       fid = fopen([fn,'_HReport.dat'],'r');

       i=1;
       while i>0
           line = fgetl(fid);
           if line==-1
               fclose(fid);
               break
           else 
               RPT{i}=strtrim(line); 
               i=i+1;
           end
       end

       %% Read in Header and Assign Parameters
       fid=fopen(f,'r');
       
       % Discard Description and first parameter line (no longer used)
       fgetl(fid); fgetl(fid);
      
       % Get heat flow pulse parameter line
       headerlines=fullfile(pwd,[fn,'_Header.dat']);
       fid1 = fopen(headerlines,'r');
       eval(['!awk ''/@ HFP/'' ',headerlines,' > ','HFParams.tmp'])
       HFParams_line=fullfile(pwd,'HFParams.tmp');
       fid2 = fopen(HFParams_line,'r');
       HFParams_line=fgetl(fid2); 
        
       % Get paramaters out
       % Remove the leading part of the line to isolate the JSON string
        % Locate the position of the JSON start
        jsonStart = strfind(HFParams_line, '{');
        jsonString = HFParams_line(jsonStart:end);
        
        % Decode the JSON string into a MATLAB structure
        parameters = jsondecode(jsonString);
        
        % Access individual parameters if needed
        loginterval = parameters.loginterval;
        heatpulseinterval = parameters.heatpulseinterval;
        power = parameters.power;
        workingdepth = parameters.workingdepth;
        resetdepth = parameters.resetdepth;
        postpenetrationdelay = parameters.postpenetrationdelay;
        accomminterval = parameters.accomminterval;
        timeafterheatpulse = parameters.timeafterheatpulse;
        autonomousmode = parameters.autonomousmode;
        accerationdeadband = parameters.accerationdeadband;
        heaterresistanceohms = parameters.heaterresistanceohms;
        autodepthmode = parameters.autodepthmode;
        headunits = parameters.headunits;
        hfp_timeafterpen = parameters.hfp_timeafterpen;

        Params = {'loginterval', 'heatpulseinterval', 'power','workingdepth',...
            'resetdepth','postpenetrationdelay','accomminterval','timeafterheatpulse',...
            'autonomousmode','accerationdeadband','heaterresistanceohms','autodepthmode',...
            'headunits','hfp_timeafterpen'};
        ParamsVals = [loginterval, heatpulseinterval, power,workingdepth,...
            resetdepth,postpenetrationdelay,accomminterval,timeafterheatpulse,...
            autonomousmode,accerationdeadband,heaterresistanceohms,autodepthmode,...
            headunits,hfp_timeafterpen];

       fclose(fid);
       fclose(fid1);
       fclose(fid2);

       % Store Header in Structure PARAMETERS
       Parameters = struct('Field',{Params},'Value',ParamsVals);

       %% Remove Temporary Files Generated by AWK
           % Delete temp file
           if exist('RawData.tmp','file')
               delete('RawData.tmp');
           end

           % Delete temp file
           if exist('RawDataTemp.tmp','file')
               delete('RawDataTemp.tmp');
           end

            % Delete temp file
           if exist('HFParams.tmp','file')
               delete('HFParams.tmp');
           end

           % Delete temp file
           if exist([fn,'_HReportTemp.dat'],'file')
               delete([fn,'_HReportTemp.dat']);
           end

       %% Move DAT files to subfolder
       movefile *.dat outputs


       
                