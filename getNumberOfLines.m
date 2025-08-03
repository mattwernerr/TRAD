function n = getNumberOfLines(fid)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Get the number of lines in a file by passing the file ID (fid) and
% counting the number of newline characters (\n) from the current pointer
% position. 
% 
% This function DESTROYS the current position of the pointer and resets it
% back to the beginning.
% 
%    Inputs:
% 
%               fid - File identifier that specifies the location in the
%                     file being read.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%    Outputs:
%                 n - Number of (remaining) lines in the file. The total
%                     number of lines in the file is obtained (ie n = N)
%                     when the file ID is passed directly after opening the
%                     file up.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 

% Set the pointer to the file's beginning
frewind(fid)

% Begin the counter
n = 0;

% Count the number of (remaining) lines the file
while (~feof(fid))
    fgetl(fid);
    n = n + 1;
end

% Reset the pointer to the file's beginning
frewind(fid)