%% Vector Field for Panagou Paper
clear all
T = 40; % simulation time
dt= 0.01; % step time
N = T/dt;
t=linspace(0,T,N);

q0 = [-10;10;pi/4]; % Arbitrary start position
v0 = [1;0.1];

Q(:,1)=q0;
V(:,1)=v0;

q = q0;
v = v0;

% System Properties

m=10;
r=0.5;
J=m*r^2/2;

sigma=1;

Qlqr = (1/2)*[1 0 0 0 0 ; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 sigma 0; 0 0 0 0 sigma];
R = (1/2)*[sigma/m^2 0; 0 sigma/J^2];

for i=2:N
    
    A = [0 0 -v(1)*sin(q(3)) cos(q(3)) 0 ;...
        0 0 v(1)*cos(q(3)) sin(q(3)) 0 ;...
        0 0 0 0 1;...
        0 0 0 0 0;....
        0 0 0 0 0];
    
    B = [0 0;
        0 0;
        0 0;
        1/m 0;
        0 1/J];
    
    [K,S,E] = lqr(A,B,Qlqr,R,0);
    
    u=-K*[q;v];
    
    
    vdot = [1/m 0 ; ...
        0 1/J]*u;
    
    v = v + vdot*dt;
    
    qdot = [cos(q(3)) 0 ; sin(q(3)) 0; 0 1]*v;
    
    q = q+dt*qdot;
    
    Q(:,i) = q;
    V(:,i) = v;
    U(:,i) = u;
    
    if rank(ctrb(A,B)) < 5
        printf('Uh Oh')
    end
    
end
%%
figure(1)
set(gcf,'color','w');

for i=N
    clf
    subplot(2,2,[2])
    plot(t(2:i),Q(:,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[4])
    plot(t(2:i),V(1,2:i),t(2:i),V(2,2:i),t(2:i),U(1,2:i),t(2:i),U(2,2:i));
    legend({'$v$','$\omega$','$u_1$','$u_2$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[1 3])
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Q(1,2:i),Q(2,2:i),'--');
    hold on
    plotAgent(Q(1,i),Q(2,i),Q(3,i),r,r/2)
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    axis equal
    pause(10^-6)
    hold off
end


