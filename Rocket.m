% Rocket
clear
clf
T = 40; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

% Agent Properties

m = 1; % Mass
l = 1; % length
g = 9.81 % Gravity

Ig = 1/12*(m*(2*l)^2); % MoI about centre of mass

q0=[10;-10;10000;0;0;0;0;0;0;0]; % Initial conditions

F = -40*ones(1,N);

PHI1 = -0.001*sin(1*t);
PHI2 = 0.001*sin(1*t);

q=q0;
Q(:,1)=q0;

for i=2:N
    
    theta1=q(4);
    theta2=q(5);
    
    phi1=PHI1(i);
    phi2=PHI2(i);
    
    f = F(i);
    qdot = [q(6:10) ; (-f*sin(phi1)*cos(theta2)-(f*sin(phi2)*sin(theta1)-f*cos(phi1)*cos(phi2)*cos(theta1))*sin(theta2))/m ; ...
        (f*sin(phi2)*cos(theta1)+f*cos(phi1)*cos(phi2)*sin(theta1))/m ;...
        (-f*sin(phi1)*sin(theta2)+(f*sin(phi2)*sin(theta1)-f*cos(phi1)*cos(phi2)*cos(theta1))*cos(theta2))/m - g;...
        ((-f*sin(phi1)*cos(theta2)-(f*sin(phi2)*sin(theta1)-f*cos(phi1)*cos(phi2)*cos(theta1))*sin(theta2))*l*cos(theta1)-(-f*sin(phi1)*sin(theta2)+(f*sin(phi2)*sin(theta1)-f*cos(phi1)*cos(phi2)*cos(theta1))*cos(theta2))*l*sin(theta1)*cos(theta2))/Ig;...
        (-(f*sin(phi2)*cos(theta1)+f*cos(phi1)*cos(phi2)*sin(theta1))*l*cos(theta2)-(-f*sin(phi1)*sin(theta2)+(f*sin(phi2)*sin(theta1)-f*cos(phi1)*cos(phi2)*cos(theta1))*cos(theta2))*l*cos(theta1)*sin(theta2))/Ig];
    
    q = q+dt*qdot;
    
        
    Q(:,i) = q;

end
%%
figure(1)
set(gcf,'color','w');


plot(t,Q(1:5,:))
legend('x','y','z','\theta_1','\theta_2')

%%
for i=1:5*dt*N:N
    clf
    subplot(2,2,[1 3])
    plot3(Q(1,1:i),Q(2,1:i),Q(3,1:i),'--',Q(1,i),Q(2,i),Q(3,i),'o');
    xlim([min(Q(1,:)) max(Q(1,:))])
    ylim([min(Q(2,:)) max(Q(2,:))])
    zlim([min(Q(3,:)) max(Q(3,:))])
    box on
    grid on
    subplot(2,2,2)
    plot(t(1:i),Q(:,1:i))
    xlim([0 T])
    ylim([min(min(Q)) max(max(Q))])
    legend({'$x_1$','$x_2$','$x_3$'},'Interpreter','latex','Location','best')
    box on
    subplot(2,2,4)
    plot(t(1:i),U1(1:i),t(1:i),U2(1:i))
    xlim([0 T])
    ylim([min(min(U1),min(U2)) max(max(U1),max(U2))])
    legend({'$u_1$','$u_2$'},'Interpreter','latex','Location','best')
    box on
    pause(0.00001)
end
%%
plot(t,H)
xlim([0 T])
ylim([-1 1])