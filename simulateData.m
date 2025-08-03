function [t, X, q, wb] = simulateData(s, Fa, Mb, m, cb, Ib, T)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Generate simulated data given the force F and torque T about a rigid
% body's center of mass or about a point p fixed on the rigid body which is
% not the center of mass.
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~                             START CASES                             ~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% =========================================================================
%                    Case 1: About the center of mass.
% =========================================================================
% 
% In this (easy) case, the position and velocity are determined by Newton's
% second law,
%                                 F = m*a,
% 
% where F is the external applied force, m is the total mass, and a is the
% inertial acceleration, while the angular velocity is determined by
% Euler's second law,
%                         M = I*wbdot + wb x (I*wb),
% 
% where M is the external applied torque, I is the moment of inertia matrix
% about the rigid body's center of mass expressed in the body-fixed frame
% (so that I is a constant 3-by-3 matrix), and wb is the angular velocity
% of the body-fixed frame relative to the inertial frame expressed in the
% body-fixed coordinate frame.
% 
% Note that Newton's law is 2nd order while Euler's law is 1st order.
%   - To keep Newton's law simple, express F in the inertial frame so that
%     the derivative relationship between position and velocity is trivial
%     and no ficticious forces involving the angular velocity appear in the
%     expression for the inertial acceleration.
%   - To keep Euler's law simple, express M in the body-fixed frame so that
%     the inertia matrix is constant in time. Doing this avoids using
%     rotation matrices furnished by the quaternion q.
% 
% =========================================================================
%     Case 2: About a point on the body that is NOT the center of mass.
% =========================================================================
% 
% This case is more difficult since Newton's and Euler's laws are now
% coupled together. This implies that both F and M must be expressed in the
% same frame, which also couples in the quaternion dynamics. Thus, the
% entire 13-dimensional state is generally coupled into each of the
% individual components. Due to the coupling, this case requires inverting
% a time-varying 6-by-6 matrix to solve for the linear and angular
% accelerations, which is more expensive to evaluate than case 1, and it
% also runs the risk of the matrix becoming singular during integration,
% though this event is (probably) unlikely. The relation between linear and
% angular acceleration in this case are given by
%    _ _     _                   _  _ _     _       _  _             _
%   |   |   |                     || . |   |         ||               |
%   | F |   |   m       -m*[c]x   || v |   |  1    0 || m*[w]x*[w]x*c |
%   |   | = |                     || . | + |         ||               |
%   | M |   | m*[c]x   (I-m*[c]x) || w |   | [c]x  1 ||   [w]x*I*w    |
%   |_ _|   |_                   _||_ _|   |_       _||_             _|
% 
% where [c]x and [w]x indicate the cross-product equivalent matrix. Note
% that both matrices are 6-by-6; as such, m is m*eye(3) in the top left and
% 1 = eye(3). Likewise, 0 = zero(3). As previously stated, vdot = dv/dt and
% wdot = dw/dt indicate the INERTIAL linear and angular accelerations,
% respectively. As such, every other quantity must also be expressed in the
% inertial frame, leverging the quaternion's rotation matrix.
% 
% Note that Newton's law is 2nd order while Euler's law is 1st order.
%   - To keep Newton's law "simple," express both F and M in the inertial
%     frame so that the derivative relationship between position and
%     velocity is trivial and no ficticious forces involving the angular
%     velocity appear in the expression for the inertial acceleration.
%   - Using the inertial frame implies that all body-fixed expressions
%     (specifically M, c, I, and wb) must be expressed in the inertial
%     frame. This is done by using the rotation matrix furnished by the
%     quaternion q such that M_inertial = R*M etc., except for the moment
%     of inertia which transforms tensorially, I_inertial = R*I*R'.
%   - After using the inertial frame to calculate the linear and angular
%     velocity derivatives, transform the angular velocity back to the
%     rotating frame to be consistent with the choice of state in case 1
%     and for easy compatability with the quaternion dynamics.
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~                              END CASES                              ~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% In both cases, the dynamics of the quaternion q = (x,y,z,w) are entirely
% kinematic and are coupled to the angular velocity expressed in the
% body-fixed frame (wb) via
%                                                      _           _
%                                                     |  w  -z   y  |
%                                                     |             |
%                                                     |  z   w  -x  |
%             qdot = 0.5 * g(q) * wb,    where g(q) = |             |
%                                                     | -y   x   w  |
%                                                     |             |
%                                                     |_-x  -y  -z _|
% 
% The initial conditions for the entire state x = (X, q, wb) are taken to
% be the origin, with the exception of the quaternion which is dynamically
% constrained to the unit 3-sphere. The quaternion's initial condition is
% the unit quaternion, (0,0,0,1). All together, the initial condition is
% 
%                              X(0) = (0,0,0,0,0,0),
%                             wb(0) = (0,0,0).
%                              q(0) = (0,0,0,1),
% 
%    Inputs:
% 
%                 s - Times at which the force F and torque M are applied.
%                     Explicitly, F = F(s) and M = M(s) are given
%                     (numerical) functions of time.
%                     Size: N-by-1 (vector)
%                     Units: s (seconds)
% 
%                Fa - Applied external force expressed in the inertial
%                     frame.
%                     Size: N-by-3 (matrix)
%                     Units: N (Newtons)
% 
%                Mb - Applied external torque expressed in the body-fixed
%                     frame.
%                     Size: N-by-3 (matrix)
%                     Units: N*m (Newton-meters)
% 
%                 m - Total CONSTANT mass of the rigid body.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%                cb - CONSTANT location of point p relative to the center
%                     of mass, expressed in the body-fixed frame.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%                Ib - Total CONSTANT moment of inertia of the rigid body
%                     about the center of mass expressed in the body-fixed
%                     frame.
%                     Size: 3-by-3 (matrix)
%                     Units: N*m2 (Newton-meters squared)
% 
%                 T - Total amount of time to integrate.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%    Outputs:
% 
%                 t - Time vector at which the solutions from solving the
%                     ODEs are INTERPOLATED. This time t simply returns the
%                     input s to ensure that the simulated data are in step
%                     with the given force and torque data.
%                     Size: N-by-1 (vector)
%                     Units: s (seconds)
% 
%                 X - Inertial state vector, X = (position, velocity), as
%                     solved through Newton's second law (possibly coupled
%                     with Euler's second law if the point of interest on
%                     the rigid body is not the center of mass).
%                     Size: N-by-3 (matrix)
%                     Units: m (meters) AND m/s (meters per second)
% 
%                 q - Quaternion describing the body-fixed frame's
%                     orientation relative to the inertial frame.
%                     Size: N-by-4 (matrix)
%                     Units: -
% 
%                wb - Angular velocity of the body-fixed frame relative to
%                     the inertial frame as expressed in the body-fixed
%                     frame.
%                     Size: N-by-3 (matrix)
%                     Units: rad/s (radians per second)
% 

% Set up ODE solver
opts = odeset('AbsTol', 1e-7, 'RelTol', 1e-9);
tspan = [0, T];
IC = [0;0;0;0;0;0;0;0;0;0;0;0;1];

% Determine which solver to use based on the distance of the point from the
% body's center of mass.
if (norm(cb) < eps)
    % The point is the center of mass.
    %
    % Find the inverse of the moment of inertia matrix at the center of
    % mass in the body-fixed frame. This will be used as a pre-computed
    % value for Euler's law.
    invIb = inv(Ib);
    sol = ode113(@odefun_at_COM, tspan, IC, opts);
else
    % The point is not the center of mass.
    sol = ode113(@odefun_at_p, tspan, IC, opts);
end

% Replicate the data that gets output from VR headset, including
% interpolating the results at the times where the force/moment are given.
t = s;
sol_interp = deval(sol, s);
X = sol_interp(1:3,:)';
q = sol_interp(10:13,:)';
wb = sol_interp(7:9,:)';

%% ODE Functions
    function f = odefun_at_COM(t, x)
        % 
        % Simulate rigid-body rotation through the center of mass. In this
        % case, the translational and rotational dynamics ARE NOT coupled.
        % 
        % State: x = (X, Y, Z, VX, VY, VZ, wbx, wby, wbz, x,  y,  z,  w)
        %             1  2  3   4   5   6   7    8    9   10  11  12  13
        %
        % Fill out the dynamics according to Newton's second law, Euler's
        % second law, and quaternion kinematics.
        f = zeros(13,1);
        
        % =================== NEWTON'S SECOND LAW ===================
        % d/dt (position) = velocity
        % d/dt (velocity) = F/m
        Fa_att = pchip(s, Fa', t);
        f(1:3) = x(4:6);
        f(4:6) = Fa_att/m;
        
        % =================== EULER'S SECOND LAW ====================
        % d/dt (wb) = inv(I)*[M - wb x (I*wb)]
        Mb_att = pchip(s, Mb', t);
        f(7:9) = invIb*(Mb_att - cross(x(7:9),Ib*x(7:9)));
        
        % ================== QUATERNION KINEMATICS ==================
        f(10:13) = quaternionKinematics(x(10:13), x(7:9));
    end

    function f = odefun_at_p(t, x)
        % 
        % Simulate rigid-body rotation through the a point p that is NOT
        % the center of mass. In this case, the translational and
        % rotational dynamics ARE coupled.
        % 
        % State: x = (X, Y, Z, VX, VY, VZ, wbx, wby, wbz, x,  y,  z,  w)
        %             1  2  3   4   5   6   7    8    9   10  11  12  13
        %
        % Fill out the dynamics according to Newton's second law, Euler's
        % second law, and quaternion kinematics.
        f = zeros(13,1);
        
        % Form the rotation matrix to transform COORDINATES from the body
        % frame to the inertial frame.
        Rab = quatRot(x(10:13)')';
        
        % Calculate the current force (inertial) and torque (body-fixed)
        % being applied to the point p.
        Fa_att = pchip(s, Fa', t);
        Mb_att = pchip(s, Mb', t);
        
        % Express everything in the inertial frame
        ca = Rab*cb;
        wa = Rab*x(7:9);
        Ia = Rab*Ib*(Rab');
        Ma_att = Rab*Mb_att;
        
        % Get cross-product equivalent matrices
        cax = getCrossProductEquivalentMatrix(ca);
        wax = getCrossProductEquivalentMatrix(wa);
        % Define M as the mass matrix and p as the ficticious force vector
        % (not to be confused with the torque M, which has been
        %  appropriately labelled with its coordinate system).
        M = [m*eye(3), -m*cax; m*cax, Ia-m*cax^2];
        p = [eye(3), zeros(3); cax, eye(3)]*[m*wax^2*ca; wax*Ia*wa];
        
        % ============= NEWTON'S AND EULER'S SECOND LAWS ============
        % Calculate both linear and angular accelerations in the inertial
        % frame. This way, the derivative relationship between position and
        % velocity stays simple.
        f(1:3) = x(4:6);
        f(4:9) = M\([Fa_att; Ma_att] - p);
        
        % Rotate the angular velocity dynamics back to the body-fixed
        % frame. Note that we calculated the inertial rate of change of the
        % angular velocity, but transforming specifically the angular
        % velocity's derivative requires only the usual rotation matrix.
        f(7:9) = (Rab')*f(7:9);
        
        % ================== QUATERNION KINEMATICS ==================
        % Use the angular velocity expressed in the body-fixed rotating
        % frame.
        f(10:13) = quaternionKinematics(x(10:13), x(7:9));
    end

    function dqdt = quaternionKinematics(q, wb)
        %
        % Calculate the quaternion kinematic dynamics given the current
        % quaternion and angular velocity of the body-fixed frame relative
        % to the inertial frame expressed in the body-fixes frame at the
        % current timestep. Use the quaternion with the scalar part last,
        % i.e. q = (x,y,z,w) where w is the scalar.
        %
        
        % Calculate dq/dt
        g = [q(4)*eye(3) + getCrossProductEquivalentMatrix(q(1:3));
             -q(1:3)'];
        dqdt = 0.5*g*wb;
    end
end