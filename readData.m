function [t, X, q] = readData(file)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% 
% Read data file (.tsv) from Unity containing...
%   1. Timestamps (t)
%   2. Position (X, Y, Z)
%   3. Quaternion (x, y, z, w)
% taken from onboard sensors of AR/VR headsets.
% 
% THE RETURNED VALUES WILL BE FILTERED TO REMOVE DUPLICATE TIMESTAMPS!
% 
% This will be done by calling unique() on the time column and removing any
% corresponding rows from the position and quaternion data.
% 
% The file can be delimited by '(' and ')' separating  time from position
% and position from quaternion.
% 
%   Ex:
%       Timestamp	Head Position (x, y, z)	Head Rotation (x, y, z, w)	
%       01/31/2023 16:23:47.314	(0, 0, 0)	(0, 0, 0, 1)
% 
%   where Head Position and Head Rotation have six digits of precision.
% 
%   Inputs:
% 
%              file - Name of file containing the data to be extracted.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%   Outputs:
% 
%                 t - Time, in seconds past start.
%                     Size: N-by-1 (vector)
%                     Units: Durations
% 
%                 X - Position.
%                     Size: N-by-3 (matrix)
%                     Units: m?
% 
%                 q - Quaternion
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
% **************** ONLY UNIQUE TIMES & DATA ARE RETURNED ****************
% 

% Open the file with read 'r' permissions.
fid = fopen(file, 'r');

% Check to make sure the file is not empty
if (feof(fid))
    error("File is empty.")
end

% Get how many number of lines there are in the file. Do this to
% preallocate memory for the returns.
N = getNumberOfLines(fid);

% Subtract 1 off the total number of lines to account for the header
N = N - 1;
% Subtract 2 off the total number of lines to account for the first two
% lines
N = N - 2;

%% Begin reading lines
% Read the first line. This should be...
% 'Timestamp	Head Position (x, y, z)	Head Rotation (x, y, z, w)	'
fgetl(fid);

% Read past the first two lines...For some reason these are always (0,0,0)
% and (0,0,0,1) followed by ~0.7 second pause before taking data.
fgetl(fid);
fgetl(fid);

% Read the fourth line. This is where the data should begin.
% Also start a counter k.
line = fgetl(fid);
k = 1;

% Get the initial time (HH:MM:SS.sss) and make this a datetime.
% Do this by splitting the line by the delimiter '('
data = strsplit(line, '(');
t0 = datetime(data{1}, 'InputFormat', 'MM/dd/yyyy HH:mm:ss.SSS');

% Preallocate memory for t, X, and q. Note that this also allocates memory
% for duplicated time values. Removing them will be the last step.
t = datetime(zeros(N,1), zeros(N,1), zeros(N,1), 'InputFormat', 'MM/dd/yyyy HH:mm:ss.SSS');
X = NaN(N,3);
q = NaN(N,4);

% Fill in the data read from the first row.
t(1) = t0;
[X(1,:), q(1,:)] = XqSplit(data);

% Begin filling in the data. The times will be converted to datetimes and
% returned relative to the initial epoch t0.
while (~feof(fid))
    line = fgetl(fid);
    k = k + 1;
    data = strsplit(line, '(');
    t(k) = datetime(data{1}, 'InputFormat', 'MM/dd/yyyy HH:mm:ss.SSS');
    [X(k,:), q(k,:)] = XqSplit(data);
end

% Eliminate all duplicate times and the associated positions/quaternions.
[t, ic, ~] = unique(t);
X = X(ic,:);
q = q(ic,:);

% Make all the times relative to the initial epoch by subtracting t0.
t = t - t0;
t.Format = 'mm:ss.SSS';

% Close the file
fclose(fid);
end

function [X, q] = XqSplit(data)
    % Remove the trailing ')' delimiter and remove any spaces.
    for k = [2 3]
        data{k} = erase(data{k}, ")");
        data{k} = erase(data{k}, " ");
    end
    
    % Now split the data again using the delimiter ',' to extract the
    % numbers.
    X = cellfun(@str2num,strsplit(data{2},','));
    q = cellfun(@str2num,strsplit(data{3},','));
    
%     Xs = strsplit(data{2},',');
%     X = str2double([string(Xs{1}), string(Xs{2}), string(Xs{3})]);
%     
%     qs = strsplit(data{3},',');
%     q = str2double([string(qs{1}), string(qs{2}), string(qs{3}), string(qs{4})]);
end
