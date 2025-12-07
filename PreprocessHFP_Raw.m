function [Parameters, dataloaded, parsedtiming, HSTART, RPT, NoTherm, NumSensUsed] = ...
         PreprocessHFP_Raw(f)
% PreprocessHFP_Raw
%
% Unified loader for all known HFP/SlugPen RAW file formats:
%   - Legacy $ format
%   - @ format with simple header
%   - @ HFP: {JSON} modern firmware
%
% Inputs:
%   f (string) - full path to .raw file
%
% Outputs:
%   Parameters     - struct of probe config parameters
%   dataloaded     - struct containing ACC, Z, PWR, TRAW, YMD/HMS
%   parsedtiming   - struct with YR MO DY HR MN SC
%   HSTART         - cell array of heat-pulse events
%   RPT            - cell array of RPT events
%   NoTherm        - number of thermistors
%   NumSensUsed    - same value (kept for SlugPen compatibility)
% Unified HFP RAW loader — supports legacy $, mid-generation @, modern JSON HFP


% -------------------------------------------------------------------------
% If no filename is given, ask the user via uigetfile()
% -------------------------------------------------------------------------
if nargin < 1 || isempty(f)
    [filename, pathname] = uigetfile({'*.raw;*.txt','Raw Probe Files'}, ...
        'Select a RAW file');
    if isequal(filename,0)
        error('No file selected.');
    end
    f = fullfile(pathname, filename);
end

% -------------------------------------------------------------------------
% From here down, you run your universal parser (format detection, etc.)
% -------------------------------------------------------------------------

%% --- Read file into memory -----------------------------------------------
fid = fopen(f,'r');
if fid < 0
    error('Cannot open file: %s', f);
end

lines = {};
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    lines{end+1,1} = line; %#ok<AGROW>
end
fclose(fid);

if isempty(lines)
    error('File is empty.');
end

%% --- Detect format (legacy $, or JSON/@) ---------------------------------
fmt = detectFormat(lines);

switch fmt
    case "json"
        Parameters = parseJsonParameters(lines);
    case "legacy"
        Parameters = parseLegacyParameters(lines);
    otherwise
        error('Unknown RAW file format.');
end

%% --- Parse the data block, RPT, HSTART -----------------------------------
[dataloaded, parsedtiming, HSTART, RPT, NoTherm] = parseData(lines);

NumSensUsed = NoTherm;

[pathstr, fname, ~] = fileparts(f);

outdir = fullfile(pathstr, 'outputs');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

matfilename = fullfile(outdir, [fname '_raw.mat']);

save(matfilename, ...
    'Parameters', ...
    'dataloaded', ...
    'parsedtiming', ...
    'HSTART', ...
    'RPT', ...
    'NoTherm', ...
    'NumSensUsed', ...
    'f', ...
    'fname');

fprintf('Saved MAT file:\n   %s\n', matfilename);

end
%% ========================================================================
%% ---------------- HELPER FUNCTIONS --------------------------------------
%% ========================================================================

function fmt = detectFormat(lines)
    fmt = "unknown";
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if isempty(s), continue; end
        if startsWith(s,"$")
            fmt = "legacy";
            return;
        elseif startsWith(s,"@")
            fmt = "json";
            return;
        end
    end
end

%% ---------------- JSON / modern format header ---------------------------
function Parameters = parseJsonParameters(lines)
% Parse JSON-style HFP header from @ HFP: {...} line
% Returns Parameters.Field (cellstr) and Parameters.Value (numeric row)

    % Find the line containing the HFP JSON
    jsonLine = '';
    for i = 1:numel(lines)
        if contains(lines{i}, '@ HFP:')
            jsonLine = lines{i};
            break;
        end
    end

    if isempty(jsonLine)
        warning('No JSON HFP header found; returning empty Parameters.');
        Parameters = struct('Field', {}, 'Value', []);
        return;
    end

    % Extract JSON substring
    idx = strfind(jsonLine, '{');
    if isempty(idx)
        warning('HFP line found but no "{" for JSON; returning empty Parameters.');
        Parameters = struct('Field', {}, 'Value', []);
        return;
    end

    jsonString = jsonLine(idx(1):end);
    p = jsondecode(jsonString);   % struct

    % ---- Explicit list of parameters we care about (numeric) ----
    Params = {'loginterval', ...
              'heatpulseinterval', ...
              'power', ...
              'workingdepth', ...
              'resetdepth', ...
              'postpenetrationdelay', ...
              'accomminterval', ...
              'timeafterheatpulse', ...
              'autonomousmode', ...
              'accerationdeadband', ...  % keep misspelling to match old code
              'heaterresistanceohms', ...
              'autodepthmode', ...
              'headunits', ...
              'hfp_timeafterpen'};

    ParamsVals = nan(size(Params));

    for k = 1:numel(Params)
        fld = Params{k};

        if isfield(p, fld)
            val = p.(fld);

        elseif strcmp(fld, 'accerationdeadband') && isfield(p, 'accelerationdeadband')
            % handle spelling difference between JSON and old code
            val = p.accelerationdeadband;

        else
            val = NaN;  % missing field → leave as NaN
        end

        % Only accept numeric/logical; otherwise NaN
        if isnumeric(val) || islogical(val)
            ParamsVals(k) = double(val);
        else
            ParamsVals(k) = NaN;
        end
    end

    Parameters = struct('Field', {Params}, 'Value', ParamsVals);
end


%% ---------------- Universal data parser ---------------------------------
function [dataloaded, parsedtiming, HSTART, RPT, NoTherm] = parseData(lines)

YMD=[]; HMS=[]; ACCx=[]; ACCy=[]; ACCz=[]; PWR=[]; Z=[]; Traw=[];
HSTART={}; RPT={};

NoTherm = [];

for i = 1:numel(lines)
    s = strtrim(lines{i});
    if isempty(s), continue; end

    % skip header lines
    if startsWith(s,'@') || startsWith(s,'$')
        continue;
    end

    % heat-pulse events
    if contains(s,'HSTART')
        HSTART{end+1} = s; %#ok<AGROW>
        continue;
    end

    % reports (skip J-)
    if contains(s,'RPT')
        if ~contains(s,'J-')
            RPT{end+1} = s; %#ok<AGROW>
        end
        continue;
    end

    % data lines require ACC and THM
    if ~contains(s,'ACC ') || ~contains(s,'THM')
        continue;
    end

    % tokenize
    tk = strsplit(s);
    if numel(tk) < 12
        continue;
    end

    % parse date
    ymd = sscanf(strrep(tk{1},'-',' '),'%d %d %d');
    if numel(ymd)~=3, continue; end

    % parse time
    hms = sscanf(strrep(tk{4},':',' '),'%d %d %d');
    if numel(hms)~=3, continue; end

    YMD(end+1,:) = ymd(:)'; %#ok<AGROW>
    HMS(end+1,:) = hms(:)'; %#ok<AGROW>

    ACCx(end+1,1) = str2double(strrep(tk{6},'X','')); %#ok<AGROW>
    ACCy(end+1,1) = str2double(strrep(tk{7},'Y','')); %#ok<AGROW>
    ACCz(end+1,1) = str2double(strrep(tk{8},'Z','')); %#ok<AGROW>
    PWR(end+1,1)  = str2double(strrep(tk{9},'V','')); %#ok<AGROW>
    Z(end+1,1)    = str2double(strrep(tk{10},'D','')); %#ok<AGROW>

    temps = str2double(tk(12:end));
    temps = temps(:);

    if isempty(NoTherm)
        NoTherm = numel(temps);
        Traw = NaN(NoTherm,0);
    end

    if numel(temps) ~= NoTherm
        temps = padOrTrim(temps, NoTherm);
    end

    Traw(:,end+1) = temps; %#ok<AGROW>
end

if isempty(NoTherm)
    error('No thermistor data found.');
end

% reorder as original SlugPen does
Traw(NoTherm+1,:) = Traw(1,:);
Traw = Traw(2:end,:);

dataloaded = struct('YMD',YMD,'HMS',HMS,'ACCx',ACCx,'ACCy',ACCy,...
                    'ACCz',ACCz,'PWR',PWR,'Z',Z,'TRAW',Traw);

parsedtiming = struct('YR',YMD(:,1),'MO',YMD(:,2),'DY',YMD(:,3), ...
                      'HR',HMS(:,1),'MN',HMS(:,2),'SC',HMS(:,3));
end

%% Pad or trim vector to length N
function x = padOrTrim(x,N)
if numel(x) > N
    x = x(1:N);
else
    x(end+1:N) = NaN;
end
end


