clear all
% Example for Lyapunov, Filipov in a Switching Scenario

T = 10; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

b=0;
k=1;
m=1;

q0 = 5*randn(2,1); % Initial Conditions (make this arbitrary)

q=q0;

Q(:,1)=q0;


for i=2:N
    qdot = [q(2);-(k/m)*q(1)-b/m*sign(q(2))];
    q = q+dt*qdot;
    Q(:,i) = q;
end

figure(2)
plot(t(1:i),Q(1,1:i),t(1:i),Q(2,1:i));
ylim([min(min(Q)),max(max(Q))]);
xlim([min(t),max(t)]);
legend({'$x~[m]$','$\dot{x}~[m/s]$'},'Interpreter','latex','Location','bestoutside')
box on


%%
i=1;
for x1 = 1.1*min(Q(1,:)):0.5:1.1*max(Q(1,:))
    for x2= 1.1*min(Q(2,:)):0.5:1.1*max(Q(2,:))
        
        dot(:,i)=[x2;-k/m*x1-b/m*sign(x2)];
        
        X(:,i)=[x1;x2];
        
        
        i=i+1;
    end
end

figure(1)
set(gcf,'color','w');

for i=1:N/100:N
    clf
    
    subplot(2,2,[1 3]);
    
    plot(Q(1,1:i),Q(2,1:i),'--r',Q(1,i),Q(2,i),'xr',[-b/m b/m],[0 0],'--k')
    hold on
    quiver(X(1,:),X(2,:),dot(1,:),dot(2,:),'k')
    
    xlim([1.1*min(Q(1,:)),1.1*max(Q(1,:))]);
    ylim([1.1*min(Q(2,:)),1.1*max(Q(2,:))]);
    xlabel('$x~[m]$','Interpreter','latex')
    ylabel('$\dot{x}~[m]$','Interpreter','latex')
    box on
    subplot(2,2,2);
    
    viscircles([Q(1,i),0],0.5,'color','k');
    hold on
    plot([1.3*min(Q(1,:))-2 Q(1,i)-0.5],[-.3 -.3],'--k','LineWidth',2)
    plot([1.3*min(Q(1,:))-2 Q(1,i)-0.5],[0.3 0.3],'-.k','LineWidth',2)
    plot([0 0],[-2 2],'--k','LineWidth',0.5)
    legend('Spring','Damper','Equilibrium','Location','best')
    ylim([-2,2]);
    xlim([1.3*min(Q(1,:))-2,1.3*max(Q(1,:))]);
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



