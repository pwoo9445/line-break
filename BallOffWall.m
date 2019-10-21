clear Q X dot
% Example for Lyapunov, Filipov in a Switching Scenario for a ball bouncing

T = 8; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

k=0.72; % Coefficient of restitution
g=9.810; % Gravitational Constant

L = 1.5; % distance to wall
h = 0; % relative height

V0 = 6.37 ; % relase velocity
theta0 = pi/3; % release angle

q0 = [0;V0*cos(theta0);h;V0*sin(theta0)]; % Initial Conditions (mak)

q=q0;

Q(:,1)=q0;

for i=2:N
    
    if abs(q(1)) <= L
        qdot = [q(2);0;q(4);-g];
        q = q+dt*qdot;
        VEL=q(2);
    elseif abs(q(1)-L)<=0.01
        q = [L;-k*q(2);q(3);q(4)];
        
        
    end
    
    Q(:,i) = q;
    
end
%%
clf
figure(2)
for i=1:5:N
    clf
set(gcf,'color','w');
plot(Q(1,1:i),Q(3,1:i));
hold on
    plot((Q(1,i)), (Q(3,i)),'.k','MarkerSize',20);

ylim([min(min(Q(3,:))),max(max(Q(3,:)))]);
xlim([0,L]);
legend({'$x~[m]$','$\dot{x}~[m/s]$'},'Interpreter','latex','Location','bestoutside')
box on
pause(0.0001)
end

%%

clear x y
clf


th = (L)/(V0*cos(theta0));
tg = (L)/(k*V0*cos(theta0));

T=th+tg;

V0 = 3%(1+k)*L/(2*sin(theta0)*cos(theta0));

N = T/dt;
t=linspace(0,T,N);

YH = q0(3)+V0*sin(theta0)*th-9.81*th^2/2;

for i=1:N
    if t(i) <= th
        x(i) = q0(1)+V0*cos(theta0)*t(i);
        y(i) = q0(3)+V0*sin(theta0)*t(i)-9.81*t(i).^2/2;

    else
        x(i) = (L-k*V0*cos(theta0)*(t(i)-th));
        y(i) = YH-k*V0*sin(theta0)*(t(i)-th)-9.81*(t(i)-th)^2/2;

    end
end


plot(x,y)


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



