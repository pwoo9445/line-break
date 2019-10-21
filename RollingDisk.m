clear Q t constraint Z
% figure(1)
clf
% figure(2)
% clf

T = 20; % simulation time
dt= 0.01; % step time
N = T/dt;
t=linspace(0,T,N);

vc=1; % 1 if using potential function, 0 if not.

q0 = [2;3;pi/6;0;0;0;0;0]; % Arbitrary start position
q0(6)=q0(5)*tan(q0(3));
Q(:,1)=q0;
q = q0;

% Final Coordinates
xf = 0;
yf = 0;
thetaf = 0;

% System Properties

m=1;
r=1;
Jspin=m*r^2/4;
Jroll=m*r^2/2;

origin=0;

epsilon=0.05;

a=1; % Damping coefficient

K0=2;

% Potential Field
Field = 1;
% 1 = Minima, 2 = Canyon, 3 = Sine Canyon,
% 4 = Gaussian at (xf,yf), 5 = Dual Parabaloid, 6 = Parallel Park
% 7 = Slanted Saddle (Gate), 8 = x, 9 = Sin-Cos Repulsion Field

for i=2:N
    
    K=[K0];
    
    if thetaf > pi
        thetaf = thetaf - 2*pi;
    end
    
    % Artificial Potential Function V(x,y,theta) and Derivative
    if Field == 1
        V = (1/2)*(q(1)^2 + q(2)^2 + q(3)^2); % Minima
        dV = [1*(q(1)-xf);1*(q(2)-yf);0*(q(3)-thetaf);0]; % Minima
    elseif Field == 2
        slope=1;
        V = -(q(1)^2 - slope*q(2) + q(3)^2); % Canyon
        dV = [2*(q(1)-xf);-slope*1;0;0]; % Canyon
    elseif Field == 3
        V = -((q(2)-sin(q(1)))^2 - q(1)); % Sine canyon
        dV = [2*cos(q(1))*(q(2)-sin(q(1)))-2;2*(q(2)-sin(q(1)));0;0]; % Sine canyon
    elseif Field == 4
        sigma = norm(q0(1)-xf,q0(2)-yf);
        dV = [(1/(2*sigma^4*pi)*(q(1)-xf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));(1/(2*sigma^4*pi)*(q(2)-yf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));0;0];
    elseif Field == 5
        n=2; % Directional quadratic. Even powers.
        m=2; % Braking polynonmial
        m1 = tan(thetaf);
        m2 = tan(thetaf-pi/2);
        if thetaf==-pi/2 || pi/2
            dV = [(n/2)*(q(1)-xf)^(n-1);(m/2)*(q(2)-yf)^(m-1);0;0];
        else
            dV = [(n/2)*m1*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*m2*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);(n/2)*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);0;0];
        end
    elseif Field == 6
        dV = -[sech(pi-1*q0(1)/2)^2;2*q(2);0;0];
    elseif Field == 7
        K=K0;
        dV = 4*[q(1)/(q(2)^2+1);(-(q(1)^2-1)*q(2)+(q(2)^4)/2+2*q(2)^2+1)/(q(2)^2+1)^2;0;0];
    elseif Field == 7
        K=K0;
        dV = 4*[q(1)/(q(2)^2+1);(-(q(1)^2-1)*q(2)+(q(2)^4)/2+2*q(2)^2+1)/(q(2)^2+1)^2;0;0];
    end
    
    % Dissipative Force
    Fd = -a*[q(5);q(6);q(7);q(8)];
    
    %     A = [1/sqrt(m+r^2*cos(q(3))^2*Jroll) 0 0 -r*cos(q(3))/sqrt(m+r^2*cos(q(3))^2*Jroll) ; 0 1/sqrt(m+r^2*sin(q(3))^2*Jroll) 0 -r*sin(q(3))/sqrt(m+r^2*sin(q(3))^2*Jroll)];
    A = [1 0 0 -r*cos(q(3)) ; 0 1 0 -r*sin(q(3))];
    D = [r*cos(q(3)) 0 ; r*sin(q(3)) 0 ; 0 1 ; 1 0  ];
    
    G = [m 0 0 0; 0 m 0 0; 0 0 Jspin 0; 0 0 0 Jroll];
    Gsharp = inv(G);
    
    % Orthogonal Projections
    
    Pd =  (D)*inv(transpose(D)*D)*transpose(D);
    Pdperp = transpose(A)*inv(A*transpose(A))*A;
    
    %     Pdperpdot=[2*r^2*cos(q(3))*(q(7))*sin(q(3))/(r^2+1), r^2*(q(7))*sin(q(3))^2/(r^2+1)-r^2*cos(q(3))^2*(q(7))/(r^2+1), 0, r*(q(7))*sin(q(3))/(r^2+1);...
    %         r^2*(q(7))*sin(q(3))^2/(r^2+1)-r^2*cos(q(3))^2*(q(7))/(r^2+1), -2*r^2*cos(q(3))*(q(7))*sin(q(3))/(r^2+1), 0, -r*(q(7))*cos(q(3))/(r^2+1);...
    %         0, 0, 0, 0;...
    %         r*(q(7))*sin(q(3))/(r^2+1), -r*(q(7))*cos(q(3))/(r^2+1), 0, 0];
    %     Pd+Pdperp;
    
    %     % Calculate Virtual Torques
        PHI = [(dV(1)/sqrt(dV(1)^2+1));(dV(2)/sqrt(dV(2)^2+1));(dV(3)/sqrt(dV(3)^2+1));(dV(4)/sqrt(dV(4)^2+1))];
        Fstar = K*PHI;
        dellPF = [0 0 -sin(q(3)) 0 ;...
            0 0 cos(q(3)) 0 ;...
            0 0 0 0 ;...
            0 0 0 0 ];...
    
        Qstar = transpose(dellPF)*Fstar;
    
    Lambda = [r*q(7)*q(8)*sin(q(3));-r*q(7)*q(8)*cos(q(3)); 0 ;0 ];
    
    qdot = [q(5) ; q(6) ; q(7) ; q(8); Pd*(-dV+Gsharp*(Fd+Qstar))-Lambda];%+[0;0;0;0;Qstar];
    
    q = q+dt*qdot;
    
    constraint(:,i) =  A*qdot(1:4); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    %     Fin = (Fd(1)+dV(1))*sin(q(3))+(Fd(2)+dV(2))*cos(q(3));
    %     tau = F(3);
    %
    Q(:,i) = q;
    U(:,i) = [-q(3)*Jspin;-m*r*(q(1)*cos(q(3))+q(2)*sin(q(3)))];
    L(:,i)=Lambda;
    %     Lstar(:,i)=Lstar;
    
    %     if abs(q(1)) < epsilon && abs(q(2)) < epsilon && abs(q(3)) < epsilon
    %         origin = origin+1;
    %     end
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

x = linspace(min(Xmin,xf)-xstep,max(Xmax,xf+xstep));
y = linspace(min(Ymin,yf-ystep),max(Ymax,yf+ystep));
[X,Y] = meshgrid(x,y);
if Field==1
    Z = 1/2*(1*(X-xf).^2+1*(Y-yf).^2); % Minima
elseif Field==2
    Z = (X.^2-slope*Y); % Canyon
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
    Z = tanh(X-1*q0(1)/2)+1+Y.^2; % G
elseif Field==7
    Z = (X.^2+Y.^2)./(1+Y.^2)+Y;
end

Zmax = max(max(abs(Z)));
Z=1*Z./Zmax;

figure(1)
set(gcf,'color','w');
clf
for i=N:500:N
%     clf
    subplot(2,2,[2])
    plot(t(2:i),Q(1:8,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$','$\phi$','$\dot{x}$','$\dot{y}$','$\dot{\theta}$','$\dot{\phi}$'},'Interpreter','latex','Location','best')
    grid on
    box on
    subplot(2,2,[4])
    plot(t(2:i),L(:,2:i),t(2:i),constraint(:,2:i),t(2:i),U(:,2:i));
    legend({'$\lambda_1$','$\lambda_2$','$\alpha_1 \dot{q}$','$\alpha_2 \dot{q}$','$\tau_1$','$\tau_2$'},'Interpreter','latex','Location','best'); %
    grid on
    box on
    
    subplot(2,2,[1 3])
    %     surf(X,Y,Z-1);
    %     shading interp;
    %     colormap jet;
    %     colorbar('Ticks',[min(min(Z))-1,0],...
    %         'TickLabels',{'Minimum','Maximum'})  ;
    %     az = 0; el = 90; view(az, el);
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    hold on
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Q(1,2:i),Q(2,2:i),'k--','Linewidth',1);
    hold on
    %     viscircles([0 0],epsilon);
    %     viscircles([0 0],norm(q0(1:2)),'Color','b');
    plotDisk(Q(1,i),Q(2,i),Q(3,i),-Q(4,i),r/2)
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    set(gca,'ZTick',[])
    grid off
    box on
    pause(10^-9)
    hold off
end



