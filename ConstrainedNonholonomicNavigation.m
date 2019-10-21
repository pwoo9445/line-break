clear Q t constraint Z PHI Lambda
% figure(1)
clf
% figure(2)
% clf

T = 20; % simulation time
dt= 0.02; % step time
N = T/dt;
t=linspace(0,T,N);

vc=1; % 1 if using potential function, 0 if not.

q0 = [1;0;pi/2;0;0;0]; % Arbitrary start position
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

G = [m 0 0 ; 0 m 0 ; 0 0 J];
Gsharp = inv(G);

% Control Variables

a=2; % Damping coefficient

K0=4;

% Navigation Properties

Nw = 11; % Number of waypoints


w = 3*randn(Nw,2); % Random coordinates of waypoints
% circular path
zeta = linspace(0,2*pi,Nw);
w = [cos(zeta)' sin(zeta)'];

w(1,:)=[q(1) q(2)]; % First waypoint is initial coordinate
W = length(w); % Number of waypoints
s=1; % Slope
rw=0.3  ; % Waypoint radius
b=1;
j=2; % Waypoint counter

% Obstacles
No = 0; % Number of obstacles
o = [3 3]; % Random coordinates of waypoints
lo=1; % obstacles width
wo=1; % obstacle length
thetao=pi/4;
Ao=0;
sigma_x=(wo*sin(thetao)+lo*cos(thetao));
sigma_y=(wo*cos(thetao)+lo*sin(thetao));
% ellipsoid properties
a = 2*wo;
b = 1*lo;

for i=2:N
    clf
    K=K0;%sign(q(1)-xf)*sign(q(2)-yf)*sign(q(3)-thetaf)*K0*tanh(norm([q(1)-xf,q(2)-yf,q(3)-thetaf])); %-sign(q(1)-xf)*sign(q(2)-yf)*
    
    xmin=min(w(j-1,1),w(j,1));
    xmax=max(w(j-1,1),w(j,1));
    
    ymin=min(w(j-1,2),w(j,2));
    ymax=max(w(j-1,2),w(j,2));
    
    x=linspace(xmin-rw,xmax+rw,100);
    y=linspace(ymin-rw,ymax+rw,100);
    
    dy = (w(j,2)-w(j-1,2));
    dx = (w(j,1)-w(j-1,1));
    
    m = dy/dx;
    
    alpha = atan(s);
    theta = atan2(dy,dx);
    
    [X,Y] = meshgrid(x,y);
    
    xo=5-0*t(i);
    yo=xo;
    
    if dx == 0
        Z = (1/(2*r))*(X-w(j-1,1)).^2-s*(X*cos(theta)+Y*sin(theta))+Ao*exp(-((X-xo).^2/(2*sigma_x^2)+(Y-yo).^2/(2*sigma_y^2)));
        dV = -[1/r*(q(1)-w(j-1,1)).^1-s*cos(theta);-s*sin(theta)];
        x = [w(j,1) w(j,1)];
        y = [w(j-1,2) w(j,2)];
        
    else % Gaussian Object
        Z = (1/b)*1*(1/(2*r))*(Y-(m*(X-w(j,1))+w(j,2))).^2-1*s*(X*cos(theta)+Y*sin(theta))+Ao*exp(-((X-xo).^2/(2*sigma_x^2)+(Y-yo).^2/(2*sigma_y^2)));
        dV = -[(1/b)*1/r*(-m)*(q(2)-(m*(q(1)-w(j,1))+w(j,2))).^1-1*s*cos(theta)-(Ao*exp(-((q(1)-xo).^2/(2*sigma_x^2)))*2*((q(1)-xo)/(2*sigma_x^2)));
            (1/b)*1/r*(q(2)-(m*(q(1)-w(j,1))+w(j,2))).^2-1*s*sin(theta)-(Ao*exp(-((q(2)-yo).^2/(2*sigma_y^2)))*2*((q(2)-yo)/(2*sigma_y^2))) ]; % Canyon
        x = xmin-r:0.1:xmax+r;
        y = m*(x-w(j,1))+w(j,2);
        %     else % Ellipsoid Obbject
        %         Z = (1/b)*1*(1/(2*r))*(Y-(m*(X-w(j,1))+w(j,2))).^2-1*s*(X*cos(theta)+Y*sin(theta))+Ao*real(sqrt(1-(((X-xo).*cos(theta)+(Y-yo).*sin(theta))/a).^2-(((Y-yo).*cos(theta)-(X-xo).*sin(theta))/b).^2));
        %         dV = -[(1/b)*1/r*(-m)*(q(2)-(m*(q(1)-w(j,1))+w(j,2))).^1-1*s*cos(theta)-Ao*((q(1)-xo)/(a^2*sqrt(1-(((q(1)-xo)*cos(theta)+(q(2)-yo)*sin(theta))/a)^2-(((q(2)-yo)*cos(theta)-(q(1)-xo)*sin(theta))/b)^2)));
        %             (1/b)*1/r*(q(2)-(m*(q(1)-w(j,1))+w(j,2))).^2-1*s*sin(theta)-Ao*((q(2)-yo)/(b^2*sqrt(1-(((q(1)-xo)*cos(theta)+(q(2)-yo)*sin(theta))/a)^2-(((q(2)-yo)*cos(theta)-(q(1)-xo)*sin(theta))/b)^2))) ]; % Canyon
        %         x = xmin-r:0.1:xmax+r;
        %         y = m*(x-w(j,1))+w(j,2);
    end
    
    
    dV = [dV;0];
    
    % System Distribution and Co-distribution;
    A = [-sin(q(3)) cos(q(3)) 0];
    D = [cos(q(3)) 0 ; sin(q(3)) 0 ; 0 1];
    
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
    
    F = Pd*(dV+Fd+Qstar);
    
    qdot = [q(4) ; q(5) ; q(6) ; Gsharp*(F)+Lambda];
    
    q = q+dt*qdot;
    
    constraint(:,i) =  Pdperp*q(4:6);%q(4)*sin(q(3))-q(5)*cos(q(3)); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    Fin = (Fd(1)+dV(1))*sin(q(3))+(Fd(2)+dV(2))*cos(q(3));
    tau = F(3);
    
    Q(:,i) = q;
    U(:,i) = [Fin;tau];
    L(:,i)=Lambda;
    
    if (norm(w(j,:).'-q(1:2))) <= rw
        j = j + 1;
    end
    
    %     Z = Z./max(max((Z)))-(j-1)*1;
    
    Q(:,i) = q;
    
    plot(x,y,w(:,1),w(:,2),'b.');
    hold on
    surface(X,Y,Z-5)
    plotAgent(Q(1,i),Q(2,i),Q(3,i),r,r/2)
    line([q(1) (q(1)+dV(1))],[q(2) (q(2)+dV(2))],'color','r')
    plot(w(1,1),w(1,2),'kX')
    plot(Q(1,2:i),Q(2,2:i),'k--','Linewidth',1);
    hold off
    axis equal
    viscircles(w,rw*ones(1,length(w)));
%     viscircles([xo yo],lo);
    
    pause(0.000001);
end

% figure(1)
% set(gcf,'color','w');
% clf
%
%
% Xmin=min(Q(1,:))-2;
% Ymin=min(Q(2,:))-2;
% Xmax=max(Q(1,:))+2;
% Ymax=max(Q(2,:))+2;
%
% xstep = Xmax - Xmin;
% ystep = Ymax - Ymin;
%
% x = linspace(xf-xstep,xf+xstep);
% y = linspace(yf-ystep,yf+ystep);
% [X,Y] = meshgrid(x,y);
% if Field==1
%     Z = 1/2*(1*(X-xf).^2+1*(Y-yf).^2); % Minima
% elseif Field==2
%     Z = (X.^2-1*Y); % Canyon
% elseif Field==3
%     Z = (Y-sin(X)).^2-2*X; % Sine Canyon
% elseif Field==4
%     Z = 1-(1/2*sigma^2*pi)*exp(-((X-xf).^2+(Y-yf).^2)/(2*sigma^2));
% elseif Field==5
%     if thetaf==-pi/2 || pi/2
%         Z = (X-xf).^n+(Y-yf).^m;
%     else
%         Z = ((Y-yf)+m1*(X-xf)).^n+((Y-yf)+m2*(X-xf)).^m; % Ellitpic Parabaloid
%     end
% elseif Field==6
%     Z = tanh(X-1*q0(1))+1+Y.^2; % G
% elseif Field==7
%     Z = (X.^2+Y.^2)./(1+Y.^2)+Y;
% end
%
% Zmax = max(max(abs(Z)));
% Z=5*Z./Zmax;
%%
for i=N:2000:N
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

%%
figure(2)
% plot(t,Lambda1-Lambda2+Lambda1-Lambda3+Lambda2-Lambda3)
clf
plot(t,Lambda1-Lambda3)



