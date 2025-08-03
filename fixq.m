function q = fixq(q)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Apply corrections to the quaternion q by finding where there are jumps or
% discontinuities (if any) in the signal and multiply everything in between
% by -1. Do this to ensure continuity of the quaternions.
% 
% This is a legal operation since, assuming Unity is still correctly
% tracking orientation, orientation is described by both q and -q.
% 
%    Inputs:
% 
%                 q - Time series of a quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%    Outputs:
% 
%                 q - Time series of a quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 

%% Checks
% Make sure there are no NaN values
if (any(isnan(q),'all'))
    error("Quaternion contains NaN values. Quitting...")
end

%% Find Jumps
% Calculate |dq| and see where the difference between neighboring points is
% large. Do this by calling diff() and then take the norm. It's not
% possible for all four quaternions to be near 0; in fact if q --> -q, then
% there will be exactly a distance of 2 between q and -q. Since we are
% stepping in time, this distance will be slightly less than two. Thus, we
% will consider any spikes greater than 1.5 as discontinuities.
normdq = vecnorm(diff(q),2,2);
normdq(normdq < 1.5) = 0;
ic = find(normdq > 1.5);
for k = 1:2:numel(ic)-1
    ick = ic(k)+1;
    ickp1 = ic(k+1);
    q(ick:ickp1,:) = -q(ick:ickp1,:);
end

% If there were an odd number of discontinuties found, then the last one
% did not get picked up by the for loop. Since the last one doesn't have a
% pair, let's assume that it trails through the very end. Thus, flip
% everything from the last discontinuity to the end.
if (k == numel(ic) - 2)
    % Ignore the possibility that k is near the very end already.
    ick = ic(k+2)+1;
    q(ick:end,:) = -q(ick:end,:);
end

% Swap
for k = 1:size(q,1)
    % Assign the current quaternion to a temporary variable and reorganize
    % the components. The scalar part always stays the same because it is
    % simply cos(Phi/2) - i.e. not attached to a coordinate system and the
    % rotation angle Phi is in an even function. Thus, only q(1:3) can
    % change. Copy the coordinate change, (3,-1,2) --> (1,2,3), which
    % changes its handedness (left-handed to right-handed), so also reverse
    % the direction of the rotation angle, Phi --> -Phi.
    qktmp  =  q(k,1:3);
    q(k,1) = -qktmp(3);
    q(k,2) =  qktmp(1);
    q(k,3) = -qktmp(2);
end

% Take q(t) relative to q(0) at time t = 0
%
% In MATLAB's quatdivide(), the scalar part goes first, so have to do some
% shifting around...
tmpq = quatdivide([q(:,4), q(:,1:3)], [q(1,4), q(1,1:3)]);
q = [tmpq(:,2:4), tmpq(:,1)];