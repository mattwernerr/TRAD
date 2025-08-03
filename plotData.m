function plotData(t,X,q)
%%
figure(1)
subplot(2,1,1)
plot(t,X)
ylabel("Position")
legend("$x$", "$y$", "$z$", 'interpreter', 'latex')
subplot(2,1,2)
plot(t,q)
ylabel("Quaternion")
legend("$x$", "$y$", "$z$", "$w$", 'interpreter', 'latex')

%%
figure(2)
scatter3(q(:,1), q(:,2), q(:,3), 1, q(:,4))
cb = colorbar;
cb.Limits = [-1,1];
cb.Ticks = -1:0.2:1;
xlabel("$x$")
ylabel("$y$")
zlabel("$z$")
ylabel(cb,"$w$", 'Interpreter', 'latex')

figure(3)
plot3(X(:,1), X(:,2), X(:,3))
%% Animation
figure(4)
R = quatRot(q);
ei = [1;0;0];
ej = [0;1;0];
ek = [0;0;1];
% % Plot a smiley face directly in front
% th = linspace(-pi,pi,50);
% cmplxCircle = exp(1i*th);
for k = 1:size(R,3)
%     subplot(2,4,[1,2,5,6])
    plot3([0,ei(1)],[0,ei(2)],[0,ei(3)],'r'), text(1.1,0,0,"$x$",'fontsize',15,'Interpreter','latex');
    hold on
    plot3([0,ej(1)],[0,ej(2)],[0,ej(3)],'g'), text(0,1.1,0,"$y$",'fontsize',15,'Interpreter','latex')
    plot3([0,ek(1)],[0,ek(2)],[0,ek(3)],'b'), text(0,0,1.1,"$z$",'fontsize',15,'Interpreter','latex')
    xlim([-3,3])
    ylim([-3,3])
    zlim([-3,3])
    grid on
    
    % Calculate new basis vector directions
    C = R(:,:,k);
    Rei = C*ei;
    Rej = C*ej;
    Rek = C*ek;
    
    % Plot the new basis, following the position coordinates as well
    Xk = X(k,:);
    xk = Xk(1);
    yk = Xk(2);
    zk = Xk(3);
    plot3(xk+[0,Rei(1)], yk+[0,Rei(2)], zk+[0,Rei(3)],'r');
    plot3(xk+[0,Rej(1)], yk+[0,Rej(2)], zk+[0,Rej(3)],'g');
    plot3(xk+[0,Rek(1)], yk+[0,Rek(2)], zk+[0,Rek(3)],'b');
    
%     % Plot a smiley face
%     plot3(3*ones(size(th)), real(0.25*cmplxCircle), 2+imag(0.25*cmplxCircle),'.--k')
%     plot3(3, +0.05, 2.05,'.b')
%     plot3(3, -0.05, 2.05,'.b')
%     kk = floor(numel(th)/4);
%     plot3(3*ones(size(th(1:kk))), real(0.125*cmplxCircle(1:kk)), 2+imag(0.125*cmplxCircle(1:kk)),'.--k')
    
    axis equal
    xlim([-3,3])
    ylim([-3,3])
    zlim([-3,3])
    hold off
    
    view(120,22)
%     view(-Rei)
    
%     subplot(2,4,[3,4])
%     plot(t,X)
%     xlabel("$t$")
%     ylabel("Position")
%     xline(t(k), 'linewidth',1)
%     legend("$x$", "$y$", "$z$", 'interpreter', 'latex', 'Location', 'eastoutside')
%     subplot(2,4,[7,8])
%     plot(t,q)
%     xlabel("$t$")
%     ylabel("Quaternion")
%     xline(t(k), 'linewidth',1)
%     legend("$x$", "$y$", "$z$", "$w$", 'interpreter', 'latex', 'Location', 'eastoutside')
    pause(0.01)
%     if (k == 1)
%         gif('90clockwisethenLeft(quaternions_normalized).gif','DelayTime',1/60);
%     elseif (mod(k,2) == 0)
%         gif;
%     end
end

figure(5)
smoothed_wb = smoothdata(wb, 1, 'sgolay', 15);
plot(t,wb/(2*pi),'.')
hold on
plot(t,smoothed_wb(:,1)/(2*pi), 'color', '#0072BD')
plot(t,smoothed_wb(:,2)/(2*pi), 'color', '#D95319')
plot(t,smoothed_wb(:,3)/(2*pi), 'color', '#EDB120')
hold off
ylabel("Angular Velocity [rev/s]")
h=legend("$\omega_x$", "$\omega_y$", "$\omega_z$", "$\omega_x$", "$\omega_y$", "$\omega_z$", 'interpreter', 'latex', 'numcolumns', 2, 'location', 'best');
title(h, "Raw $\mid$  Filtered")