clear Q t constraint Z PHI Lambda Edot
% figure(1)
clf
% figure(2)
% clf

T = 20; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

vc=1; % 1 if using potential function, 0 if not.

q0 = [-5;5;0;0;0;0]; % Arbitrary start position
q0(5)=q0(4)*tan(q0(3));
Q(:,1)=q0;
q = q0;

% Final Coordinates
xf = 0;
yf = 0;
thetaf = 0;

% System Properties

m=1;
r=1;
J=m*r^2/2;

origin=0;

epsilon=0.05;

a=2; % Damping coefficient

K0=10;

% Potential Fields
Field = 3;
% 1 = Minima, 2 = Canyon, 3 = Sine Canyon,
% 4 = Gaussian at (xf,yf), 5 = Dual Parabaloid, 6 = Parallel Park
% 7 = Slanted Saddle (Gate)

for i=2:N
    
    K=K0;%sign(q(1)-xf)*sign(q(2)-yf)*sign(q(3)-thetaf)*K0*tanh(norm([q(1)-xf,q(2)-yf,q(3)-thetaf])); %-sign(q(1)-xf)*sign(q(2)-yf)*
    %     K = K0*tanh(thetaf-q(3));
    if thetaf > pi
        thetaf = thetaf - 2*pi;
    end
    
    % Artificial Potential Function V(x,y,theta) and Derivative
    if Field == 1
        V = -(1/2)*(q(1)^2 + q(2)^2); % Minima
        dV = -[1*(q(1)-xf);1*(q(2)-yf);(0*q(3)-thetaf)]; % Minima
        
    elseif Field == 2
        V = -(q(1)^2 - 1*q(2) + q(3)^2); % Canyon
        dV = -[2*(q(1)-0);-1*1;0]; % Canyon
    elseif Field == 3
        V = -((q(2)-sin(q(1)))^2 - q(1)); % Sine canyon
        dV = -[2*cos(q(1))*(q(2)-sin(q(1)))-2;2*(q(2)-sin(q(1)));0*q(3)]; % Sine canyon
        K=K0;
    elseif Field == 4
        sigma = norm(q0(1)-xf,q0(2)-yf);
        dV = -[(1/(2*sigma^4*pi)*(q(1)-xf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));(1/(2*sigma^4*pi)*(q(2)-yf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));0];
    elseif Field == 5
        n=2; % Directional quadratic. Even powers.
        m=2; % Braking polynonmial
        m1 = tan(thetaf);
        m2 = tan(thetaf-pi/2);
        if thetaf==-pi/2 || pi/2
            dV = -[(n/2)*(q(1)-xf)^(n-1);(m/2)*(q(2)-yf)^(m-1);0];
        else
            dV = -[(n/2)*m1*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*m2*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);(n/2)*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);0];
        end
    elseif Field == 6
        K0=2;
        a=1;
        K=-sign(q(1)-xf)*sign(q(2)-yf)*K0;
        dV = -[sech(pi-1*q0(1))^2;2*q(2);0];
    elseif Field == 7
        K=K0;
        dV = -[q(1)/(q(2)^2+1);(-(q(1)^2-1)*q(2)+(q(2)^4)/2+2*q(2)^2+1)/(q(2)^2+1)^2;0];
    end
    
    % System Distribution and Co-distribution;
    A = [-sin(q(3)) cos(q(3)) 0];
    D = [cos(q(3)) 0 ; sin(q(3)) 0 ; 0 1];
    Adot = [-cos(q(3)) -sin(q(3)) 0];
    
    G = [m 0 0 ; 0 m 0 ; 0 0 J];
    Gsharp = inv(G);
    
    Pd =  (D)*((transpose(D)*D)\transpose(D));
    Pdperp = transpose(A)*((A*transpose(A))\A);
    
    
    Pdperpdot = [2*cos(q(3))*(q(6))*sin(q(3)), (-2*cos(q(3))^2+1)*(q(6)), 0;...
        (-2*cos(q(3))^2+1)*(q(6)), -2*cos(q(3))*(q(6))*sin(q(3)), 0;...
        0, 0, 0];
    
    % Dissipative Force
    Fd = -a*[q(4);q(5);q(6)];
    
    % Virtual Surface
    PHI = [(dV(1)/sqrt(dV(1)^2+1));(dV(2)/sqrt(dV(2)^2+1));(dV(3)/sqrt(dV(3)^2+1))];
    
    % Virtual Torque
    Fstar = K*PHI;
    dellPF = [0 0 -sin(q(3)) ;...
        0 0 cos(q(3))  ;...
        0 0 0 ];...
        
    Qstar = transpose(dellPF)*Fstar;
    % Base lambda
    Lambda = -Pdperpdot*q(4:6); % From dynamic equations
    
    
    %     Lambda1(:,i) = -Pdperpdot*q(4:6); % From dynamic equations
    %     Lambda2(:,i) = inv(A*transpose(A))*transpose(A)*q(6)*(q(4)*cos(q(3))+q(5)*sin(q(3)));
    %     Lambda3(:,i) = Adot*q(4:6); % From dynamic equations
    
    %     % Virtual Lambda
    %     Lstar = L - inv(A*Gsharp*transpose(A))*A*Gsharp*(dV+Fd);
    %
    %     F = dV + Fd + transpose(A)*(Lstar-L);
    
    % Control Forces (Orthogonal Projections of surface and dissipative
    % force)
    
    F = Pd*(dV+Fd);
    
    qdot = [q(4) ; q(5) ; q(6) ; Gsharp*(F+Qstar)+Lambda];
    
    q = q+dt*qdot;
    
    Edot(i) = transpose(Gsharp*(F+Qstar)+Lambda)*q(4:6); % Lyapunov Stability. Check if Edot = F * gamma' < 0 for all q.
    
    constraint(:,i) =  Pdperp*q(4:6);%q(4)*sin(q(3))-q(5)*cos(q(3)); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    Fin = (Fd(1)+dV(1))*cos(q(3))+(Fd(2)+dV(2))*sin(q(3));
    tau = F(3);
    
    Q(:,i) = q;
    U(:,i) = [Fin;tau];
    L(:,i)=Lambda;
    
    if abs(q(1)) < epsilon && abs(q(2)) < epsilon && abs(q(3)) < epsilon
        origin = origin+1;
    end
end

figure(1)
set(gcf,'color','w');
clf


Xmin=min(Q(1,:))-2;
Ymin=min(Q(2,:))-2;
Xmax=max(Q(1,:))+2;
Ymax=max(Q(2,:))+2;

xstep = Xmax - Xmin;
ystep = Ymax - Ymin;

x = linspace(xf-xstep,xf+xstep);
y = linspace(yf-ystep,yf+ystep);
[X,Y] = meshgrid(x,y);
if Field==1
    Z = 1/2*(1*(X-xf).^2+1*(Y-yf).^2); % Minima
elseif Field==2
    Z = (X.^2-1*Y); % Canyon
elseif Field==3
    Z = (Y-sin(X)).^2-2*X; % Sine Canyon
elseif Field==4
    Z = 1-(1/2*sigma^2*pi)*exp(-((X-xf).^2+(Y-yf).^2)/(2*sigma^2));
elseif Field==5
    if thetaf==-pi/2 || pi/2
        Z = (X-xf).^n+(Y-yf).^m;
    else
        Z = ((Y-yf)+m1*(X-xf)).^n+((Y-yf)+m2*(X-xf)).^m; % Ellitpic Parabaloid
    end
elseif Field==6
    Z = tanh(X-1*q0(1))+1+Y.^2; % G
elseif Field==7
    Z = (X.^2+Y.^2)./(1+Y.^2)+Y;
end

Zmax = max(max(abs(Z)));
Z=5*Z./Zmax;

for i=1:200:N
    %     clf
    subplot(2,2,2)
    plot(t(2:i),Q(1:6,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$','$\dot{x}$','$\dot{y}$','$\dot{\theta}$'},'Interpreter','latex','Location','best')
    grid on
    box on
    subplot(2,2,4)
    plot(t(2:i),L(:,2:i),t(2:i),constraint(:,2:i),t(2:i),U(:,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    legend({'$\lambda_1$','$\lambda_2$','$\lambda_2$','${A(q) \dot{q}}_1$','${A(q) \dot{q}}_2$','${A(q) \dot{q}}_3$','$F$','$\tau$'},'Interpreter','latex','Location','best')
    grid on
    box on
    
    subplot(2,2,[1 3])
    surf(X,Y,Z-5);
    shading interp;
    colormap jet;
    colorbar('Ticks',[min(min(Z))-5,0],...
        'TickLabels',{'Minimum','Maximum'})  ;
    az = 0; el = 90; view(az, el);
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
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
    grid off
    box on
    pause(10^-9)
    hold off
end

figure(2)
set(gcf,'color','w');

plot(t,Edot)
xlabel({'$t~[s]$'},'Interpreter','latex')
ylabel({'$\dot{E} = \left<F;\gamma \prime \right> $'},'Interpreter','latex')
grid on

%%
figure(2)
% plot(t,Lambda1-Lambda2+Lambda1-Lambda3+Lambda2-Lambda3)
clf
plot(t,Lambda1-Lambda3)





