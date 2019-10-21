clear Q X dot
% Example of a pendulum on a cart.

T = 20; % simulation time
dt= 0.05; % step time
N = T/dt;
t=linspace(0,T,N);

l=2; % Length of pendulum
m=1; % Mass of Pendulum
M = 1; % Mass of cart
g=9.81; % Gravitational Constant
J=1/12*m*l^2;
q0 = [1;-1;pi;0];%pi*randn(2,1)]; % Initial Conditions (make this arbitrary)

q=q0;

Q(:,1)=q0;

poles = [-1 -2 -3 -4];

for i=2:N
    
    % Linear (Solve for control inputs)
    % This ignored where the pendulum is and just does error control
    % knowing that origin is desired trajectory
    A = [0 1 0 0; ...
        0 0 -m^2*g*l/(M*m*l^2) 0; ...
        0 0 0 1; ...
        0 0 m*g*l*(M+m)/(M*m*l^2) 0];
    
    B = [ 0 ; m*l^2/(M*m*l^2) ; 0 ; -m*l/(M*m*l^2)];
    
    K = place(A,B,poles);
    
    u = -K*q;
    
    umax =100;
    
    if u > umax
        u=umax;
    elseif u < -umax
        u=-umax;
    end
    
    % Non-linear
    %     p=(-cos(q(3))^2*m+M+m);
    %     qdot = [q(2);(m*l*q(4)^2*sin(q(3))-cos(q(3))*m*g*sin(q(3))+u)/p;q(4);((-m*l*q(4)^2*sin(q(3))-u)*cos(q(3))+(M+m)*g*sin(q(3)))/(l*p)];
    % Linear Feedback
    
    qdot = A*q + B*u
%     qdot = (A-B*K)*q;
    
    % Discrete-integral
    q = q+dt*qdot;
    
    % Add some disturbances
    
    if mod(i,round(N/6))==0
    q(1) = q(1) + 0.1*randn;
    end
    if mod(i,round(N/32))==0
    q(3) = q(3) + 0.02*randn;
    end

    % Save states and inputs 
    Q(:,i) = q;
    U(i)=u;
    
end

figure(2)
set(gcf,'color','w');

plot(t,Q);
ylim([min(min(Q)),max(max(Q))]);
xlim([min(t),max(t)]);
legend({'$x~[m]$','$\dot{x}~[m/s]$','$\theta~[rad]$','$\dot{\theta}~[rad/s]$'},'Interpreter','latex','Location','bestoutside')
box on
grid on

%%
figure(1)
set(gcf,'color','w');

for i=1:N/200:N
    clf
    subplot(2,2,[1 2]);
    plot([Q(1,i) Q(1,i)+l*sin(Q(3,i))],[0 l*cos(Q(3,i))],'k','LineWidth',2)
    hold on
    filledCircle([Q(1,i)+l*sin(Q(3,i)) l*cos(Q(3,i))],0.1,100,'k');
    rectangle('Position',[Q(1,i)-1.5 -1 3 1],'Curvature',0.2);
    filledCircle([Q(1,i)-1 -1.1],0.1,100,'k');
    filledCircle([Q(1,i)+1 -1.1],0.1,100,'k');
    plot([Q(1,i) Q(1,i)-0.05*U(i)],[-0.5 -0.5],'--r','LineWidth',1);
    plot([min(Q(1,:))-l max(Q(1,:))+l],[-1.2 -1.2],'--k','LineWidth',0.5);
    ylim([-1.5*l 1.5*l]);
    xlim([min(Q(1,:))-l max(Q(1,:))+l]);
    ylabel('$y~[m]$','Interpreter','latex')
    xlabel('$x~[m]$','Interpreter','latex')
    axis equal
    box on
    subplot(2,2,[3 4]);
    plot(t(1:i),Q(:,1:i),t(1:i),U(1:i),'--r','LineWidth',1);
    ylim([min(min(Q)),max(max(Q))]);
    xlim([min(t),max(t)]);
    legend({'$x~[m]$','$\dot{x}~[m/s]$','$\theta~[rad]$','$\dot{\theta}~[rad/s]$','u'},'Interpreter','latex','Location','best')
    box on
    pause(0.000001)
    
end



