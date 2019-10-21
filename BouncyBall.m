clear Q X dot
% Example for Lyapunov, Filipov in a Switching Scenario for a ball bouncing

T = 15; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

k=0.9; % Coefficient of restitution
g=9.81; % Gravitational Constant

h = 10; % Starting height

q0 = [h;2]; % Initial Conditions (mak)

q=q0;

Q(:,1)=q0;

for i=2:N
    
    if q(1) > 0.0
        qdot = [q(2);-g];
        q = q+dt*qdot;
        VEL=q(2);
    elseif q(1) <=0.000001
        q = [0.000001;-k*VEL];
    end
    
    Q(:,i) = q;
    
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
for x1 = 1.1*min(Q(1,:)):(1.1*max(Q(1,:))-1.1*min(Q(1,:)))/100:1.1*max(Q(1,:))
    for x2= 1.1*min(Q(2,:)):(1.1*max(Q(2,:))-1.1*min(Q(2,:)))/100:1.1*max(Q(2,:))
        if x1 > 0.0
            dot(:,j) = [x2;-g];
            VEL=q(2);
        elseif x1 <=0.0000001
            dot(:,j)  = [-k*VEL;0];
        end
        X(:,j)=[x1;x2];
        j=j+1;
    end
end

figure(1)
set(gcf,'color','w');

for i=1:N/200:N
    clf
    subplot(2,2,[1 3]);
    plot(Q(1,1:i),Q(2,1:i),'--r',Q(1,i),Q(2,i),'xr','MarkerSize',10)
    hold on
    quiver(X(1,:),X(2,:),dot(1,:),dot(2,:),'k')
    xlim([1.1*min(Q(1,:)) 1.1*max(Q(1,:))]);
    ylim([1.1*min(Q(2,:)) 1.1*max(Q(2,:))]);
    xlabel('$x~[m]$','Interpreter','latex')
    ylabel('$\dot{x}~[m]$','Interpreter','latex')
    box on
    
    subplot(2,2,2);
    plot([-2 2],[-.4 -.4],'k','LineWidth',2)
    hold on
    plot(0, (Q(1,i)),'.k','MarkerSize',20);
    ylim([-1 1.5*h]);
    ylabel('$x~[m]$','Interpreter','latex')
    box on
    
    subplot(2,2,4);
    plot(t(1:i),Q(1,1:i),t(1:i),Q(2,1:i));
    ylim([min(min(Q)),max(max(Q))]);
    xlim([min(t),max(t)]);
    legend({'$x~[m]$','$\dot{x}~[m/s]$'},'Interpreter','latex','Location','best')
    box on
    
    pause(0.000001)
    
end



