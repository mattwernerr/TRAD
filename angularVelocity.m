function wb = angularVelocity(q,qdot)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Calculate the body angular velocity according to the quaternion, q, and
% its time rate of change, qdot. Given q = (x,y,z,w), where w is the scalar
% part, do this using the kinematic relation
%                                        _ . _
%                  _                 _  |  x  |
%                 |  w   z   -y   -x  | |  .  |
%                 |                   | |  y  |
%          wb = 2 | -z   w    x   -y  | |  .  |
%                 |                   | |  z  |
%                 |_ y  -x    w   -z _| |  .  |
%                 \___________________/ |_ w _|
%                          A(q)         \_____/
%                                        dq/dt
% 
% where A(q) is a time-varying 3-by-4 matrix and dq/dt is qdot.
% 
%    Inputs:
% 
%                 q - Time-series of a quaternion.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%              qdot - Time-series of a quaternion's time derivative.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%    Outputs:
% 
%                wb - Time-series of the body-fixed angular velocity.
%                     Size: N-by-3 (matrix)
%                     Units: 1/s
% 

%% Checks
% Make sure q and qdot are the same size. This could become a problem if
% the derivative is taken using diff() or otherwise something weird happens
% along the way that removes elements.
assert(size(q,1) == size(qdot,1), "q and qdot are not the same size.");

% Fill in the angular velocity.
wb = NaN(size(q,1),3);
for k = 1:size(q,1)
    % Distribute out elements for clarity
    x = q(k,1);
    y = q(k,2);
    z = q(k,3);
    w = q(k,4);
    xdot = qdot(k,1);
    ydot = qdot(k,2);
    zdot = qdot(k,3);
    wdot = qdot(k,4);
    % Write out the kinematics equation component-wise so that the matrix
    % multiplication is already performed.
    wb(k,1) = 2*( w*xdot + z*ydot - y*zdot - x*wdot);
    wb(k,2) = 2*(-z*xdot + w*ydot + x*zdot - y*wdot);
    wb(k,3) = 2*( y*xdot - x*ydot + w*zdot - z*wdot);
end