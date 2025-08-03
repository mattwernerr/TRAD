function X = fixX(X)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Apply a coordinate rotation to the spatial variables X = (x,y,z). 
% This rotation will be LEFT-HANDED to make the results make sense in a
% more convential right-coordinate coordinate system. The rotation to be
% applied is
%                    _           _
%                   |  0   0   1  |
%                   |             |
%               A = | -1   0   0  |.
%                   |             |
%                   |_ 0   1   0 _|
% 
% This rearrangement makes X FORWARDS, Y LEFT, and Z UP, which is opposed
% to Unity's default (X RIGHT, Y UP, and Z FORWARDS)
% 
%    Inputs:
% 
%                 X - Time series of position.
%                     Size: N-by-3 (matrix)
%                     Units: -
% 
%    Outputs:
% 
%                 X - Time series of position.
%                     Size: N-by-3 (matrix)
%                     Units: -
% 

% Swap
for k = 1:size(X,1)
    % Assign the current row to a temporary variable and reorganize the
    % components so that x is forwards, y is left, and z is up.
    Xktmp  =  X(k,:);
    X(k,1) =  Xktmp(3);
    X(k,2) = -Xktmp(1);
    X(k,3) =  Xktmp(2);
end

% Take X(t) relative to X(0) at time t = 0
X = X - X(1,:);