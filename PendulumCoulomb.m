clear Q X dot
% Example for Lyapunov, Filipov in a Switching Scenario

T = 20; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

b=0;
k=1; % Friction Coefficient
l=2; % Length of pendulum
m=1; % Mass of Pendulum
g=9.81; % Gravitational Constant

q0 = pi*randn(2,1); % Initial Conditions (make this arbitrary)

q=q0;

q(1)=pi+.01;
q(2)=0;

Q(:,1)=q0;


for i=2:N
    qdot = [q(2);-(g/l)*sin(q(1))-(k/m)*(q(2))];
    q = q+dt*qdot;
    Q(:,i) = q;
    
    V(i)= 0.5*q(1)^2+0.5*q(2)^2;
    Vdot(i)= q(1) + q(2);
end

figure(2)
set(gcf,'color','w');

plot(t(1:i),Q(1,1:i),t(1:i),Q(2,1:i));
ylim([min(min(Q)),max(max(Q))]);
xlim([min(t),max(t)]);
legend({'$x~[m]$','$\dot{x}~[m/s]$'},'Interpreter','latex','Location','bestoutside')
box on


%%
j=1;
for x1 = 1.1*min(Q(1,:)):pi/10:1.1*max(Q(1,:))
    for x2= 1.1*min(Q(2,:)):pi/10:1.1*max(Q(2,:))
        
        dot(:,j)=[x2;-(g/l)*sin(x1)-(k/m)*(x2)];
        
        X(:,j)=[x1;x2];
        
%         v = 0.5*x1^2+0.5*x2^2;
%         
%       V(:,j) = [x1;x2;v];
        
        j=j+1;
    end
end
%%
figure(1)
set(gcf,'color','w');

for i=1:N/200:N
    clf
    
    subplot(2,2,[1 3]);
    plot(Q(1,1:i),Q(2,1:i),'--r',Q(1,i),Q(2,i),'xr')
    hold on
    quiver(X(1,:),X(2,:),dot(1,:),dot(2,:),'k')
%     plot3(V)
    xlim([1.1*min(Q(1,:)) 1.1*max(Q(1,:))]);
    ylim([1.1*min(Q(2,:)) 1.1*max(Q(2,:))]);
    xlabel('$\theta~[m]$','Interpreter','latex')
    ylabel('$\dot{\theta}~[m]$','Interpreter','latex')
    box on
    subplot(2,2,2);
    
    plot([0 l*sin(Q(1,i))],[0 -l*cos(Q(1,i))],'k','LineWidth',2)
    hold on
    filledCircle([l*sin(Q(1,i)) -l*cos(Q(1,i))],0.1,100,'k');
    
    ylim([-1.5*l 1.5*l]);
    xlim([-1.5*l 1.5*l]);
    ylabel('$y~[m]$','Interpreter','latex')
    xlabel('$x~[m]$','Interpreter','latex')
    box on
    subplot(2,2,4);
    plot(t(1:i),Q(1,1:i),t(1:i),Q(2,1:i));
    
    ylim([min(min(Q)),max(max(Q))]);
    xlim([min(t),max(t)]);
    legend({'$x~[m]$','$\dot{x}~[m/s]$'},'Interpreter','latex','Location','best')
    box on
    pause(0.000001)
    
end



