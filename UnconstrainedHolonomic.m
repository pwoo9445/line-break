%% Vector Field for Panagou Paper
clear Q t constraint
clf
% figure(1)
% clf
% figure(2)
% clf

T = 60; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);
 
vc=1; % 1 if using potential function, 0 if not.

q0 = [100;100;0;0;0;0]; % Arbitrary start position
q0(5)=q0(4)*tan(q0(3));
Q(:,1)=q0;
q = q0;

% System Properties

m=1;
r=1;
J=m*r^2/2;

k1=0.5;
k2=0;

origin=0;

epsilon=0.05;

alpha=m;

for i=2:N
        
    F = k1*sqrt(q(4)^2+q(5)^2);
    tau = k2*q(6);
    
    qdot = [q(4) ; q(5) ; q(6) ; 1/m*(-q(1)+F*cos(q(3))) ;1/m*(-q(2)+F*sin(q(3))) ; 1/J*(-q(3))+tau];
    
    %         figure(2)
    %         plot(t(i),q(4)*sin(q(3))+q(5)*cos(q(3)),'.')
    %         hold on
    
    q = q+dt*qdot;
    
    constraint(i) =  q(4)*sin(q(3))-q(5)*cos(q(3)); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    Q(:,i) = q;
    %     U(:,i) = u;
    L(i)=lambda;
    
    if abs(q(1)) < epsilon && abs(q(2)) < epsilon && abs(q(3)) < epsilon
        origin = origin+1;
    end
end
origin
figure(1)
set(gcf,'color','w');
clf

for i=N:10000:N
    clf
    subplot(2,2,[2])
    plot(t(2:i),Q(1:6,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$','$\dot{x}$','$\dot{y}$','$\dot{\theta}$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[4])
    plot(t(2:i),L(2:i));
    legend({'$\lambda$','$\dot{p}_1$','$p_1$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[1 3])
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Q(1,2:i),Q(2,2:i),'--');
    hold on
    viscircles([0 0],epsilon);
    viscircles([0 0],norm(q0(1:2)),'Color','b');
    
    plotAgent(Q(1,i),Q(2,i),Q(3,i),r,r/2)
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    axis equal
    pause(10^-6)
    hold off
end


