%%% ====================================================================
%%  Purpose: 
%       This function collects the input file used to load raw data 
%       in to SlugPen.
%  
%%  Last edit:
%       07/19/2023 by Kristin Dickerson, UCSC
%%% =====================================================================

function    [f, NoTherm, ...
             NumSensUsed,...
             filename] ...
             = GetFiles(figure_Main, ...
              H, ...
              edit_NumberofTherm)

%% Initialize
            % Open GUI to select filename
            % --------------------------
                [filename, pathname] = uigetfile({'*.raw';'*.mat'}, ...
                    'Pick a raw data probe output file');
                figure(figure_Main)
        
            % Validate the file/path
            % -----------------------
                if isequal(filename,0)||isequal(pathname,0)
                    uialert(figure_Main, 'File not found or not valid',...
                        'Error finding file');
                    H.Error=1;
                    return
                end

            % Get file name
            % ------------
            f=fullfile(pathname,filename);
        
            % Get the first line and confirm default header start
            fid = fopen(f,'r');

            % If fails, alert user 
            if fid<0
                uialert(figure_Main, 'Unable to open selected file','Error');
                H.Error=1;
                return

%% Load in file
            % Else, load in file    
            else
                line = fgetl(fid);
                if isempty(line)||~strcmp(line(1),'@')
                    fclose(fid);
                    disp('Invalid File Format');
                    H.Error=1;
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
                    H.Fileinfo.No_Thermistors.Value = nThermistors;  % Set the No_Thermistors in H.Fileinfo Structure
                    NoTherm = nThermistors;
                    NumSensUsed = NoTherm;
                    edit_NumberofTherm.Value = num2str(nThermistors);
                end
            end
            fclose(fid);

