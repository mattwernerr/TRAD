function [F, Mp] = ForceNTorque(M, COM_p, IMOI_COM, Xddot_p, wb, wbdot, q)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Calculate the force F and torque M at the point p, which is generally
% taken to be NOT the center of mass, given...
%   1. The characteristics of the rigid-body, including:
%       a. Mass (M)
%       b. Center of mass relative to the point p (COM_p), and 
%       c. Moment of inertia matrix about the center of mass (IMOI_COM),
% 
%   2. and point p's trajectory in the form of:
%       d. Linear acceleration (Xddot_p),
%       e. Angular velocity in the body-fixed frame (wb),
%       f. Angular acceleration in the body-fixed frame (wbdot),
%       g. Orientation (q)
% 
%    Inputs:
% 
%                 M - Total mass of the rigid body.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%             COM_p - Position of the rigid body's center of mass in the
%                     body-fixed frame relative to the point p. This POINT
%                     is constant in time, as required by rigid-body
%                     dynamics.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%          IMOI_COM - Moment of inertia matrix at the rigid body's center
%                     of mass expressed in the body-fixed frame. This MATRIX
%                     is constant in time, as required by rigid-body
%                     dynamics.
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m^2 (kilograms times square meters)
% 
%           Xddot_p - Time-series of the second time derivative of position 
%                     (linear acceleration) at the point p. This is
%                     expressed in the inertial lab frame.
%                     Size: N-by-3 (matrix)
%                     Units: m (meters)
% 
%                wb - Time-series of the body-fixed angular velocity.
%                     Size: N-by-3 (matrix)
%                     Units: 1/s (radians per second)
% 
%             wbdot - Time-series of the first time derivative of the
%                     body-fixed angular velocity (angular acceleration).
%                     Size: N-by-3 (matrix)
%                     Units: 1/s2 (radians per squared sec)
% 
%                 q - Time-series of the quaternion describing the
%                     orientation of the body-fixed frame relative to the
%                     inertial frame.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%    Outputs:
% 
%                 F - Calculated inertial external force at the point p,
%                     expressed in the body-fixed frame.
%                     Size: N-by-3 (matrix)
%                     Units: N (Newtons)
% 
%                Mp - Calculated inertial external torque at the point p,
%                     expressed in the body-fixed frame.
%                     Size: N-by-3 (matrix)
%                     Units: N*m (Newton-meters)
% 

%% Checks
% Check that M, COM, and IMOI are all constants of the right size.
if ~(numel(M) == 1 && numel(COM_p) == 3 && numel(IMOI_COM) == 9)
    error("Mass, center of mass, and inertia matrix must be constants.")
end
% Check that Xddot_p, wb, wbdot are all the same size.
if ~(all(size(Xddot_p) == size(wb)) && all(size(wb) == size(wbdot)))
    error("Time-series X, wb, and wbdot must all be the same size.")
end

%% Computations
% Allocate space for the returns
F = NaN(size(Xddot_p));
Mp = NaN(size(F));

% Calculate the force and moment about the point p. Note that everything
% needs to be expressed in the INERTIAL FRAME, including:
%   CENTER OF MASS,
%   MOMENT OF INERTIA,
%   POSITION,
%   ANGULAR VELOCITY, AND 
%   ANGULAR ACCELERATION!
% For this, we need the quaternion q to calculate the rotation matrices.
% Note that the rotation matrix to transfer between coordinate systems is
% the TRANSPOSE of the rotation matrix used to describe vectors rotating in
% the same coordinate system.
Rab = permute(quatRot(q), [2 1 3]);
for k = 1:size(F,1)
    % Get this rotation matrix to transform COORDINATES from the body frame
    % to the inertial frame.
    Rabk = Rab(:,:,k);
    % Get the center of mass/trajectories in column vector format and
    % transform them into the inertial frame.
    COM_pk = Rabk*COM_p;
    Xddot_pk = Xddot_p(k,:)'; % Already given in inertial coordinates (???)
    wak = Rabk*wb(k,:)';
    % Apply the special case of transport theorem, where the vector being
    % differentiated is the angular velocity itself. Thus, its derivative
    % does not need extra terms to properly calculate it in the inertial
    % frame. Simply take the derivative and apply the rotation matrix to
    % find its derivative (angular acceleration) in the inertial frame.
    wadotk = Rabk*wbdot(k,:)';
    % Express IMOI_COM in the inertial frame
    IMOI_COMak = Rabk*IMOI_COM*(Rabk');
    % Form the cross product equivalent matrices
    COM_pkx = getCrossProductEquivalentMatrix(COM_pk);
    wakx = getCrossProductEquivalentMatrix(wak);
    
    % Calculate F
    Fk = M*(Xddot_pk - COM_pkx*wadotk + wakx*wakx*COM_pk);
    % Calculate M at p in the inertial frame
    Mpk = M*COM_pkx*Xddot_pk + (IMOI_COMak - M*COM_pkx*COM_pkx)*wadotk + wakx*IMOI_COMak*wak + M*COM_pkx*wakx*wakx*COM_pk;
    
    % The inertial force F and inertial torque M have been calculated and
    % expressed in the inertial frame. Leave the force expressed in the
    % inertial frame, but rotate the torque to be expressed in the
    % body-fixed rotating frame.
    F(k,:)  = Fk;
    Mp(k,:) = (Rabk')*Mpk;
end