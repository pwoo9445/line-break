clear Q t constraint Z
% figure(1)
clf
% figure(2)
% clf

T = 60; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

vc=1; % 1 if using potential function, 0 if not.

q0 = [-10;10;0;0;0;0;0;0]; % Arbitrary start position
q0(6)=q0(5)*tan(q0(3));
Q(:,1)=q0;
q = q0;

% System Properties

m=1;
r=1;
I=m*r^2/2;
J=m*r^2/2;

origin=0;

epsilon=0.05;

a=2; % Damping coefficient

K=10;

Method = 1; % 1 for Control Forces, 2 for Virtual Field, Dissipative forces

% Potential Field
Field = 3; % 1 = Minima, 2 = Canyon, 3 = Sine Canyon

for i=2:N
    
    % Artificial Potential Function V(x,y,theta) and Derivative
    if Field == 1
        V = -(1/2)*(q(1)^2 + q(2)^2); % Minima
        dV = -[1*(q(1)-0);1*(q(2)-0);(0*q(3)-0)]; % Minima
        
    elseif Field == 2
        V = -(q(1)^2 - 1*q(2) + q(3)^2); % Canyon
        dV = -[2*(q(1)-0);-1*1;0]; % Canyon
        
    elseif Field == 3
        V = -((q(2)-sin(q(1)))^2 - q(1)); % Sine canyon
        dV = -[2*cos(q(1))*(q(2)-sin(q(1)))-1;2*(q(2)-sin(q(1)));0*q(3)]; % Sine canyon
    elseif Field == 4
        V = -(q(1)+q(2))*(sin(q(1)*q(2))*cos(q(1)*q(2)));
        dV = -[q(2)* (q(1) + q(2))*cos(q(1)*q(2))^2 + cos(q(1)*q(2))*sin(q(1)*q(2)) - q(2)*(q(1) + q(2))*sin((q(1)*q(2)))^2;q(1)*(q(1) + q(2))*cos(q(1)*q(2))^2+ cos(q(1)*q(2))*sin(q(1)*q(2)) - q(1)*(q(1) + q(2))*sin(q(1)*q(2))^2;0];
    end
    
    % Dissipative Force
    Fd = -a*[q(5);q(6);q(7);q(8)];
    
    A = [1 0 0 -r*cos(q(3)) ; 0 1 0 -r*sin(q(3))];
    dA = [0 0 0 r*sin(q(3)) ; 0 1 0 -r*cos(q(3))];
    G = [m 0 0 0 ; 0 m 0 0; 0 0 J 0;0 0 0 I];
    Gsharp = inv(G);
    
    % Base lambda
    L = m*q(6)*(q(4)*cos(q(3))+q(5)*sin(q(3)));
    L = inv(A*Gsharp*transpose(A))*dA*[q(5);q(6);q(7);q(8)];
    % Virtual Lambda
    Lstar = L - inv(A*Gsharp*transpose(A))*A*Gsharp*(dV+Fd);
    
    F = dV + Fd + transpose(A)*(Lstar-L);
    
    if Method == 1 % Control Forces
        qdot = [q(5) ; q(6) ; q(7) ;q(8); 1/m*(F(1)+L(1)) ; 1/m*(F(2)+L(2)) ; 1/J*(F(3)+K*(-dV(1)*sin(q(3))/sqrt(dV(1)^2+1)+dV(2)*cos(q(3))/sqrt(dV(2)^2+1)))];
    end
    
    if Method == 2 % Virtual Field and dissipative forces
        qdot = [q(4) ; q(5) ; q(6) ; 1/m*(-Lstar*sin(q(3))+dV(1)+Fd(1)) ;1/m*(Lstar*cos(q(3))+dV(2)+Fd(2)) ; 1/J*(F(3)-K*(dV(1)/sqrt(((dV(1))^2+1))+dV(2)/sqrt(((dV(2))^2+1))))];
    end
    
    q = q+dt*qdot;
    
    constraint(i) =  q(4)*sin(q(3))-q(5)*cos(q(3)); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    if mod(q(3),pi)==0
        Fin = F(1) / cos(q(3));
    else
        Fin = F(2) / sin(q(3));
    end
    tau = F(3);
    
    Q(:,i) = q;
    U(:,i) = [Fin;tau];
    L(i)=L;
    
    if abs(q(1)) < epsilon && abs(q(2)) < epsilon && abs(q(3)) < epsilon
        origin = origin+1;
    end
end
%%
origin
figure(1)
set(gcf,'color','w');
clf


Xmin=min(Q(1,:))-2;
Ymin=min(Q(2,:))-2;
Xmax=max(Q(1,:))+2;
Ymax=max(Q(2,:))+2;

x = linspace(Xmin,Xmax);
y = linspace(Ymin,Ymax);
[X,Y] = meshgrid(x,y);
if Field==1
    Z = 1/2*(1*X.^2+1*Y.^2); % Minima
elseif Field==2
    Z = (X.^2-1*Y); % Canyon
elseif Field==3
    Z = (Y-sin(X)).^2-X; % Sine Canyon
elseif Field==4
    Z = (X+Y).*sin(X.*Y).*cos(X.*Y);
end

Zmax = max(max(abs(Z)));
Z=5*Z./Zmax

for i=N:1000:N
    %     clf
    subplot(2,2,[2])
    plot(t(2:i),Q(1:6,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$','$\dot{x}$','$\dot{y}$','$\dot{\theta}$'},'Interpreter','latex','Location','best')
    grid on
    box on
    subplot(2,2,[4])
    plot(t(2:i),L(2:i),t(2:i),constraint(2:i),t(2:i),U(:,2:i));
    legend({'$\lambda$','$A(q) \dot{q}$','$F$','$\tau$'},'Interpreter','latex','Location','best')
    grid on
    box on
    
    subplot(2,2,[1 3])
    surf(X,Y,Z-5);
    shading interp;
    colormap jet;
    colorbar('Ticks',[min(min(Z))-5,0],...
        'TickLabels',{'Minimum','Maximum'})  ;
    az = 0; el = 90; view(az, el);
    xlim([Xmin Xmax]);
    ylim([Ymin Ymax]);
    hold on
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Q(1,2:i),Q(2,2:i),'k--','Linewidth',1);
    hold on
    %     viscircles([0 0],epsilon);
    %     viscircles([0 0],norm(q0(1:2)),'Color','b');
    plotAgent(Q(1,i),Q(2,i),Q(3,i),r,r/2)
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    set(gca,'ZTick',[])
    axis equal
    grid off
    box on
    pause(10^-9)
    hold off
end


