%% Simulate data
% Time from 0 to 35 seconds in 30 millisecond increments
dt = 0.030;
s = (0:dt:35)';
% Assume that the force, moment, and body characterization (mass, center of
% mass, and inertia matrix) are perfectly known (i.e., reality)
Fa = 2.5*exp(-(s-4).^2).*cos(s).*[sin(s), (s-5).^2, 0.24*sqrt(s)] + exp(-0.75*(s-10).^2).*[sin(s), cos(s), cos(2*sqrt(s))] + exp(-0.25*(s-20).^2).*[cos(s)+cos(2*s),sin(2*s)-sin(3*s),cos(3*s)+sin(2*s)];
Mb = 5*[sin(s)./(s+1), cos(s)./(s+1), erfc(s/6)] + exp(-(s-15).^2).*[sin(s+1), cos(s-2), cos(5*sin(s))] - exp(-(s-27).^2).*[0.5*s,-0.2*s.^0.6,log(s+1)]/10;
Fa = Fa;
Mb = Mb;
m = 5;
cb = [-0.1;0;-0.75];
% c = zeros(3,1);
Ib = [1.052027097046708, 0.028102471692480, 0.040637141530706
      0.028102471692480, 1.066161164789907, 0.011437850854736
      0.040637141530706, 0.011437850854736, 1.147720858123187]*10;

% Simulate the output response from the VR headset
[t_sim, X_sim, q_sim, wb_sim] = simulateData(s, Fa, Mb, m, cb, Ib, max(s));
% Take derivatives of the 'data' and calculate the angular velocity
% (expressed in the body frame)
for K = 1:1
[Xdot_sim, qdot_sim, Xddot_sim, qddot_sim] = XqDerivatives(X_sim,q_sim,t_sim,'sgolay', 15);
wb_rec = angularVelocity(q_sim,qdot_sim);
% Calculate the angular acceleration (expressed in the body frame)
for k = 3:-1:1
    wbdot_rec(:,k) = gradient(wb_rec(:,k), t_sim);
end

% Introduce errors from reality, induced by taking measurements.
e_m = 0.001*0;
e_c = 0.05*0;
e_I = 0.05*0;
error1x1 = 1 + (2*rand(1,1)-1)*e_m;
error3x1 = 1 + (2*rand(3,1)-1)*e_c;
error3x3 = 1 + (2*rand(3,3)-1)*e_I;

% Based off of the 'data' alone, calculate our guess for what we think the
% original force & moment were during the VR 'test.'
[F_rec, M_rec] = ForceNTorque(error1x1.*m, error3x1.*cb, error3x3.*Ib, Xddot_sim, wb_rec, wbdot_rec, q_sim);

subplot(1,2,1)
plot(s, Fa, '.') % Plot the exact force
hold on
% Plot the guess of the exact force
plot(t_sim, F_rec(:,1),'color','#0072BD')
plot(t_sim, F_rec(:,2),'color','#D95319')
plot(t_sim, F_rec(:,3),'color','#EDB120')
hold off
xlabel("$t$")
ylabel("Force")
h=legend("$F_X$", "$F_Y$", "$F_Z$", "$F_X$", "$F_Y$", "$F_Z$", 'interpreter', 'latex', 'numcolumns', 2, 'location', 'northeast');
title(h, "Given $\mid$  Calculated")
ylim([-1.6,1.6])

subplot(1,2,2)
plot(s, Mb, '.') % Plot the exact torque
hold on
% Plot the guess of the exact torque
plot(t_sim, M_rec(:,1),'color','#0072BD')
plot(t_sim, M_rec(:,2),'color','#D95319')
plot(t_sim, M_rec(:,3),'color','#EDB120')
hold off
xlabel("$t$")
ylabel("Torque")
h=legend("$M_x$", "$M_y$", "$M_z$", "$M_x$", "$M_y$", "$M_z$", 'interpreter', 'latex', 'numcolumns', 2, 'location', 'northeast');
title(h, "Given $\mid$  Calculated")
ylim([-1.6,1.6])
pause(0.01)
% if (K == 1)
%     gif('simulatedData_reconstruction_ForceTorque_5_percent_error_view2.gif', 'DelayTime',1)
% else
%     gif
% end
end