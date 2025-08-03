function xdot = takeDerivative(x, t)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Take a derivative of the N-dimensional function x(t), where 
% x : R --> RN.
% 
%   Inputs:
% 
%             t - Total time in seconds.
%                 Size: N-by-1 (array of Durations or doubles)
%                 Units: s
% 
%             x - Time-history of N-dimensional function, where each
%                 component of x goes across the rows and time flows down
%                 each column.
%                 Size: N-by-M (matrix)
%                 Units: ?
% 
%    Outputs:
% 
%          xdot - Time-derivative of x.
%                 Size: N-by-M (matrix)
%                 Units: ?/s
% 


% Calculate the total amount of seconds from the given amount of duration.
if (isa(t, 'duration'))
    s = seconds(t);
elseif (isa(t, 'double'))
    s = t;
else
    error("Time t must be an array of durations or doubles.")
end

% Take the derivative by calling gradient() and don't do any further
% processing. Both x(:,k) and t are N-by-1 vectors, so there is no
% ambiguity which direction gradient() is operating along.
for k = size(x,2):-1:1
    xdot(:,k) = gradient(x(:,k), s);
end
