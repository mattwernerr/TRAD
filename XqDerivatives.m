function [Xdot, qdot, Xddot, qddot] = XqDerivatives(X,q,t,filter,window)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Take the first and second derivatives of the the position X and
% quaternion q in time.
% 
%    Inputs:
% 
%                 X - Time-history of position.
%                     Size: N-by-3 (matrix)
%                     Units: ?
% 
%                 q - Time-history of quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%                 t - Time from the initial epoch.
%                     Size: N-by-1 (array of Durations) OR array of doubles
%                     Units: s
% 
%               **************** OPTIONAL ****************
% 
%               See smoothdata() doc for more information.
% 
%            filter - String providing which filter method is requested for
%                     smoothing out the derivatives.
%                     Size: 1-by-1 (string scalar)
%                     Units: N/A
%                       OPTIONS: movmean
%                                movmedian
%                                gaussian
%                                lowess
%                                loess
%                                rlowess
%                                rloess
%                                sgolay
% 
%            window - Window length, specified as a positive integer
%                     scalar, a two-element vector of positive integers, or
%                     A TWO-ELEMENT VECTOR OF POSITIVE DURATIONS.
%                     Size: 1-by-1 (scalar) OR 1-by-2 (vector)
%                     Units: -
% 
%    Outputs:
% 
%              Xdot - First derivative of position.
%                     Size: N-by-3 (matrix)
%                     Units: ?/s
% 
%              qdot - First derivative of quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: 1/s
% 
%             Xddot - Second derivative of position.
%                     Size: N-by-3 (matrix)
%                     Units: ?/s2
% 
%             qddot - Second derivative of quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: ?/s2
% 

%% Check
% Check that the number of inputs is between 3 and 5
narginchk(3, 5)
% Get the number of inputs to check if a filter and/or filter+window have
% been specified.
FILTER = false;
WINDOW = false;
switch nargin
    case 4
        FILTER = true;
    case 5
        FILTER = true;
        WINDOW = true;
end

%% Computations
% Take the first derivative of the position & quaternion. In any case, do
% this without filtering. These will be the reported Xdot and qdot.
Xdot = takeDerivative(X, t);
qdot = takeDerivative(q, t);

if (FILTER)
    % Smooth the first derivative to make the second derivative more well
    % behaved.
    if (WINDOW)
        % Do this using the specified filter AND window.
        smoothed_Xdot = smoothdata(Xdot, 1, filter, window);
        smoothed_qdot = smoothdata(qdot, 1, filter, window);
    else
        % Do this using the specified filter.
        smoothed_Xdot = smoothdata(Xdot, 1, filter);
        smoothed_qdot = smoothdata(qdot, 1, filter);
    end
    % Take the derivative of the smoothed first derivatives of X and q.
    Xddot = takeDerivative(smoothed_Xdot, t);
    qddot = takeDerivative(smoothed_qdot, t);
else
    % Take two derivatives of the raw data without using a filter.
    Xddot = takeDerivative(Xdot, t);
    qddot = takeDerivative(qdot, t);
end