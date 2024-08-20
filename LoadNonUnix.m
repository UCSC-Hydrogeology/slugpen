%%% =======================================================================
%%  Purpose:
%     This function loads data on Non-Unix systems for SlugPen.
%
%%  Last edit:
%      07/19/2023 - Kristin Dickerson, UCSC
%%% =======================================================================

function [Parameters,...
          dataloaded,...
          parsedtiming]...
          = LoadNonUnix(f, ...
                        figure_Main)

%% Initialize
          % Break apart input file 'f'
          [~,fn,ext] = fileparts(f);
          filename   = [fn ext];


           % Initialize variables
           % --------------------

           % Time Series Variables
           yr    = [];
           mo    = [];
           dy    = [];
           hr    = [];
           mn    = [];
           sc    = [];
           accx  = [];
           accy  = [];
           accz  = [];
           pwr   = [];
           z     = [];
           T0    = [];
           T1    = [];
           T2    = [];
           T3    = [];
           T4    = [];
           T5    = [];
           T6    = [];
           T7    = [];
           T8    = [];
           T9    = [];
           T10   = [];
           T11   = [];
           T12   = [];
           T13   = [];

                
           % Ensure raw file is correct format. Inform user if it is not
           % ------------------------------------------------------------
           fid  = fopen(f,'r');
           line = fgetl(fid);
           if length(line)~=21
               fclose(fid);
               uialert(figure_Main, 'This RAW file has the incorrect format.','Error');
               return
           end  

%% Read in header

           % Discard description and first parameter line (no longer used)
           % -------------------------------------------------------------
           frewind(fid);
           fgetl(fid); fgetl(fid);  

           % Get remaining header lines
           % -------------------------
           headerlines = 11;
           Params=cell(1,headerlines);
           ParamsVals=zeros(1,headerlines);  

           % Read in header
           % --------------
           for i=1:headerlines
               line = strtrim(fgetl(fid));
               a    = find(isspace(line)==1);
               Params{i}     = strtrim(line(1:a(end)));
               ParamsVals(i) = str2num(line(a(end):end)); 
           end  

           % Store Header in Structure PARAMETERS
           % --------------------------------------
           Parameters = struct('Field',{Params},'Value',ParamsVals);  

           % Inform user that data is loading
           % --------------------------------
           foo = strrep(filename,'_','-');
           h_wait = waitbar(0,['Loading ',foo]);
           set(h_wait,'name','Please Wait')
           k_wait = 1;  

%% Read in all data

           i=1; % Valid Data Counter

           while i>0
               line = fgetl(fid);
               if line == -1
                   fclose(fid);
                   break
               else
                   % Check for Data Line or Comment Line
                   a = find(isspace(line)==1);
                   if ~contains(line,'!') && length(a)==25
                       
                       % Data Line
                       v1  = line(1:a(1)-1);
                       v2  = line(a(1)+1:a(2)-1); %#ok<*NASGU>
                       v3  = line(a(2)+1:a(3)-1);
                       v4  = line(a(3)+1:a(4)-1);
                       v5  = line(a(4)+1:a(5)-1);
                       v6  = line(a(5)+1:a(6)-1);
                       v7  = line(a(6)+1:a(7)-1);
                       v8  = line(a(7)+1:a(8)-1);
                       v9  = line(a(8)+1:a(9)-1);
                       v10 = line(a(9)+1:a(10)-1);
                       v11 = line(a(10)+1:a(11)-1);
                       v12 = line(a(11)+1:a(12)-1);
                       v13 = line(a(12)+1:a(13)-1);
                       v14 = line(a(13)+1:a(14)-1);
                       v15 = line(a(14)+1:a(15)-1);
                       v16 = line(a(15)+1:a(16)-1);
                       v17 = line(a(16)+1:a(17)-1);
                       v18 = line(a(17)+1:a(18)-1);
                       v19 = line(a(18)+1:a(19)-1);
                       v20 = line(a(19)+1:a(20)-1);
                       v21 = line(a(20)+1:a(21)-1);
                       v22 = line(a(21)+1:a(22)-1);
                       v23 = line(a(22)+1:a(23)-1);
                       v24 = line(a(23)+1:a(24)-1);
                       v25 = line(a(24)+1:end);  
                       
                       % Assign
                       %     Yr Mo Dy
                       v  = str2num(strrep(v1,'-',' ')); %#ok<*ST2NM>
                       yr(i,:) = v(1); %#ok<*AGROW>
                       mo(i,:) = v(2);
                       dy(i,:) = v(3);

                       %     Hr Min Sec
                       v  = str2num(strrep(v4,':',' ')); 
                       (strrep(v4,':',' '));
                       hr(i,:) = v(1);
                       mn(i,:) = v(2);
                       sc(i,:) = v(3);

                       %     Acceleration
                       accx(i,:) = str2num(strrep(v6,'X',''));
                       accy(i,:) = str2num(strrep(v7,'Y',''));
                       accz(i,:) = str2num(strrep(v8,'Z',''));

                       %     Voltage
                       pwr(i,:)  = str2num(strrep(v9,'V',''));

                       %     Depth
                       z(i,:)    = str2num(strrep(v10,'D',''));

                       %     Thermistors
                       T0(i,:)   = str2num(v12); 
                       T1(i,:)   = str2num(v13);
                       T2(i,:)   = str2num(v14);
                       T3(i,:)   = str2num(v15);
                       T4(i,:)   = str2num(v16);
                       T5(i,:)   = str2num(v17);
                       T6(i,:)   = str2num(v18);
                       T7(i,:)   = str2num(v19);
                       T8(i,:)   = str2num(v20);
                       T9(i,:)   = str2num(v21);
                       T10(i,:)  = str2num(v22);
                       T11(i,:)  = str2num(v23);
                       T12(i,:)  = str2num(v24);
                       T13(i,:)  = str2num(v25);  
                       Traw  = [T1'; T2'; T3'; T4'; T5'; T6'; T7'; T8'; T9'; T10'; T11'; T12'; T13'; T0'];  
                      
                       % Increment
                       i=i+1;

                   end
               end

               % Update waitbar
               waitbar(k_wait/5000,h_wait);
               if k_wait==5000
                   k_wait=0;
               end
               k_wait=k_wait+1;
           end  

%% Save data in structures for access

           dataloaded = struct('YMD',ymd,'HMS', hms,'ACCx', accx,...
           ACCy',accy,'ACCz',accz,'PWR',pwr,'Z',z, 'TRAW',Traw);  
           parsedtiming = struct('YR',yr,'MO', mo,'DY', dy,...
           HR',hr,'MN',mn,'SC',sc);  
           close(h_wait);