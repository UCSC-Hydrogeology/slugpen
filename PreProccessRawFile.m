%%% =======================================================================
%%  Purpose:
%       This function preprocesses raw data file from a heat flow probe to be loaded into SlugPen. This includes:
% 
%       (a) collecting the input file used to load raw data 
%       in to SlugPen.
%       (b) loading data on Unix systems for SlugPen.
%
%%  Last edit:
%       09/04/2024 - Kristin Dickerson, UCSC
%%% =======================================================================

%% Get input files
% ---------------
% ---------------

    % Initialize
            % Open GUI to select filename
            % --------------------------
                [filename, pathname] = uigetfile({'*.raw';'*.mat';'*.txt'}, ...
                    'Pick a raw data probe output file');
        
            % Validate the file/path
            % -----------------------
                if isequal(filename,0)||isequal(pathname,0)
                    disp('ERROR: File not found or not valid')
                    return
                end

            % Get file name
            % ------------
            f=fullfile(pathname,filename);
        
            % Get the first line and confirm default header start
            fid = fopen(f,'r');

            % If fails, alert user 
            if fid<0
                disp('ERROR: Unable to open selected file');
                return

    % Load in file
            % Else, load in file    
            else
                line = fgetl(fid);
                if isempty(line)||~strcmp(line(1),'@')
                    fclose(fid);
                    disp('ERROR: Invalid File Format');
                    return
                else
                    % Read in header
                    % ---------------
                    for n=1:21
                        fgetl(fid);
                    end
            
                    % Determining NUMBER OF THERMISTORS
                    % ---------------------------------
                    % Read second line and determine number of thermistors 
                    % There are 11 columns of timing, acceleration, etc.
                    % Add 1 for bottom water sensor
                    line         = strtrim(fgetl(fid)); 
                    ncols        = length(find(isspace(line)));
                    nThermistors = ncols-11+1;
                    NoTherm = nThermistors;
                    NumSensUsed = NoTherm;
                end
            end
            fclose(fid);


if isunix
%% Load data in using UNIX
% ---------------
% ---------------

    %% Initialize
     
     % Get file name
     [~,fn,~]=fileparts(f);
        
    %% Scan using UNIX awk
                    
           
           
           % Grab all header lines
           eval(['!awk ''/@/'' ',f,' > ',fn,'_Header.dat'])
           % Grab all lines with Data
           eval(['!awk ''!/^@/'' ',f,' > ',fn,'_RawDataTemp.dat'])
           rtemp = fullfile(pwd,[fn,'_RawDataTemp.dat']);
           fid = fopen(rtemp,'r');
           eval(['!awk ''length > 180'' ',rtemp,' > ',fn,'_RawData.dat'])

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
           fid = fopen([fn,'_RawData.dat'],'r');
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
               HSTART{i}=strtrim(line); 
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
       eval(['!awk ''/@ HFP/'' ',headerlines,' > ','HFParams.dat'])
       HFParams_line=fullfile(pwd,'HFParams.dat');
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


       %% Move DAT files to subfolder
       movefile *.dat outputs

       % make a .mat file for loading into SlugPen
       rawfilename = filename(1:end-4);
       matfilename = ['outputs/',rawfilename,'_raw.mat'];
       save(matfilename, 'Parameters', ...
       'dataloaded',...
       'parsedtiming',...
       'HSTART',...
       'RPT', 'filename', 'f', 'NoTherm', 'NumSensUsed')

end

%% UPDATE NEEDED--> THE NON-UNIX CODE DOES NOT WORK YET WITH NEW RAW FILE FORMAT 09/04/2024

%else
%
%%% Load data in using NON-UNIX
%% -------------------------------
%% -------------------------------
%
%function [Parameters,...
%          dataloaded,...
%          parsedtiming]...
%          = LoadNonUnix(f, ...
%                        figure_Main)
%
%%% Initialize
%          % Break apart input file 'f'
%          [~,fn,ext] = fileparts(f);
%          filename   = [fn ext];
%
%
%           % Initialize variables
%           % --------------------
%
%           % Time Series Variables
%           yr    = [];
%           mo    = [];
%           dy    = [];
%           hr    = [];
%           mn    = [];
%           sc    = [];
%           accx  = [];
%           accy  = [];
%           accz  = [];
%           pwr   = [];
%           z     = [];
%           T0    = [];
%           T1    = [];
%           T2    = [];
%           T3    = [];
%           T4    = [];
%           T5    = [];
%           T6    = [];
%           T7    = [];
%           T8    = [];
%           T9    = [];
%           T10   = [];
%           T11   = [];
%           T12   = [];
%           T13   = [];
%
%                
%           % Ensure raw file is correct format. Inform user if it is not
%           % ------------------------------------------------------------
%           fid  = fopen(f,'r');
%           line = fgetl(fid);
%           if length(line)~=21
%               fclose(fid);
%               uialert(figure_Main, 'This RAW file has the incorrect format.','Error');
%               return
%           end  
%
%%% Read in header
%
%           % Discard description and first parameter line (no longer used)
%           % -------------------------------------------------------------
%           frewind(fid);
%           fgetl(fid); fgetl(fid);  
%
%           % Get remaining header lines
%           % -------------------------
%           headerlines = 11;
%           Params=cell(1,headerlines);
%           ParamsVals=zeros(1,headerlines);  
%
%           % Read in header
%           % --------------
%           for i=1:headerlines
%               line = strtrim(fgetl(fid));
%               a    = find(isspace(line)==1);
%               Params{i}     = strtrim(line(1:a(end)));
%               ParamsVals(i) = str2num(line(a(end):end)); 
%           end  
%
%           % Store Header in Structure PARAMETERS
%           % --------------------------------------
%           Parameters = struct('Field',{Params},'Value',ParamsVals);  
%
%           % Inform user that data is loading
%           % --------------------------------
%           foo = strrep(filename,'_','-');
%           h_wait = waitbar(0,['Loading ',foo]);
%           set(h_wait,'name','Please Wait')
%           k_wait = 1;  
%
%%% Read in all data
%
%           i=1; % Valid Data Counter
%
%           while i>0
%               line = fgetl(fid);
%               if line == -1
%                   fclose(fid);
%                   break
%               else
%                   % Check for Data Line or Comment Line
%                   a = find(isspace(line)==1);
%                   if ~contains(line,'!') && length(a)==25
%                       
%                       % Data Line
%                       v1  = line(1:a(1)-1);
%                       v2  = line(a(1)+1:a(2)-1); %#ok<*NASGU>
%                       v3  = line(a(2)+1:a(3)-1);
%                       v4  = line(a(3)+1:a(4)-1);
%                       v5  = line(a(4)+1:a(5)-1);
%                       v6  = line(a(5)+1:a(6)-1);
%                       v7  = line(a(6)+1:a(7)-1);
%                       v8  = line(a(7)+1:a(8)-1);
%                       v9  = line(a(8)+1:a(9)-1);
%                       v10 = line(a(9)+1:a(10)-1);
%                       v11 = line(a(10)+1:a(11)-1);
%                       v12 = line(a(11)+1:a(12)-1);
%                       v13 = line(a(12)+1:a(13)-1);
%                       v14 = line(a(13)+1:a(14)-1);
%                       v15 = line(a(14)+1:a(15)-1);
%                       v16 = line(a(15)+1:a(16)-1);
%                       v17 = line(a(16)+1:a(17)-1);
%                       v18 = line(a(17)+1:a(18)-1);
%                       v19 = line(a(18)+1:a(19)-1);
%                       v20 = line(a(19)+1:a(20)-1);
%                       v21 = line(a(20)+1:a(21)-1);
%                       v22 = line(a(21)+1:a(22)-1);
%                       v23 = line(a(22)+1:a(23)-1);
%                       v24 = line(a(23)+1:a(24)-1);
%                       v25 = line(a(24)+1:end);  
%                       
%                       % Assign
%                       %     Yr Mo Dy
%                       v  = str2num(strrep(v1,'-',' ')); %#ok<*ST2NM>
%                       yr(i,:) = v(1); %#ok<*AGROW>
%                       mo(i,:) = v(2);
%                       dy(i,:) = v(3);
%
%                       %     Hr Min Sec
%                       v  = str2num(strrep(v4,':',' ')); 
%                       (strrep(v4,':',' '));
%                       hr(i,:) = v(1);
%                       mn(i,:) = v(2);
%                       sc(i,:) = v(3);
%
%                       %     Acceleration
%                       accx(i,:) = str2num(strrep(v6,'X',''));
%                       accy(i,:) = str2num(strrep(v7,'Y',''));
%                       accz(i,:) = str2num(strrep(v8,'Z',''));
%
%                       %     Voltage
%                       pwr(i,:)  = str2num(strrep(v9,'V',''));
%
%                       %     Depth
%                       z(i,:)    = str2num(strrep(v10,'D',''));
%
%                       %     Thermistors
%                       T0(i,:)   = str2num(v12); 
%                       T1(i,:)   = str2num(v13);
%                       T2(i,:)   = str2num(v14);
%                       T3(i,:)   = str2num(v15);
%                       T4(i,:)   = str2num(v16);
%                       T5(i,:)   = str2num(v17);
%                       T6(i,:)   = str2num(v18);
%                       T7(i,:)   = str2num(v19);
%                       T8(i,:)   = str2num(v20);
%                       T9(i,:)   = str2num(v21);
%                       T10(i,:)  = str2num(v22);
%                       T11(i,:)  = str2num(v23);
%                       T12(i,:)  = str2num(v24);
%                       T13(i,:)  = str2num(v25);  
%                       Traw  = [T1'; T2'; T3'; T4'; T5'; T6'; T7'; T8'; T9'; T10'; T11'; T12'; T13'; T0'];  
%                      
%                       % Increment
%                       i=i+1;
%
%                   end
%               end
%
%               % Update waitbar
%               waitbar(k_wait/5000,h_wait);
%               if k_wait==5000
%                   k_wait=0;
%               end
%               k_wait=k_wait+1;
%           end  
%
%%% Save data in structures for access
%
%           dataloaded = struct('YMD',ymd,'HMS', hms,'ACCx', accx,...
%           ACCy',accy,'ACCz',accz,'PWR',pwr,'Z',z, 'TRAW',Traw);  
%           parsedtiming = struct('YR',yr,'MO', mo,'DY', dy,...
%           HR',hr,'MN',mn,'SC',sc);  
%           close(h_wait);
%
%end

