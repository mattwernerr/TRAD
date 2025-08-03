function R = quatRot(q)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Given a time series of quaternions (x, y, z, w), calculate the rotation
% matrix.
% 
%    Inputs:
% 
%                 q - Time series of quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%    Outputs:
% 
%                 R - Time series of the rotation matrix described by q.
%                     Size: 3-by-3-by-N (3D array)
%                     Units: -
% 

% Allocate space for R.
N = size(q,1);
R = NaN(3,3,N);

% Iterate through each row in q to fill out the rotation matrices one step
% at a time.
for n = 1:N
    % Split out the scalar and vector components. Also make the vector
    % component a column vector.
    qs = q(n,4);
    qv(:,1) = q(n,1:3);
    
    % Calculate the rotation matrix from this quaternion.
    qvx = getCrossProductEquivalentMatrix(qv);
    R(:,:,n) = ((qs^2 - qv'*qv)*eye(3) + 2*qv*(qv') + 2*qs*qvx);
end