% Car, Front Wheel Drive Model

clear all;

% figure(1)
clf
% figure(2)
% clf

T = 10; % simulation time
dt= 0.00001; % step time
N = T/dt;
t=linspace(0,T,N);

vc=1; % 1 if using potential function, 0 if not.

q0 = [5;5;pi/10;0;0;0;0;0]; % Arbitrary start position

Q(:,1)=q0;
q = q0;

% Final Coordinates
xf = 0;
yf = 0;
thetaf = 0;

% System Properties

m=1;    % Vehicle Mass
l=2;    % vehicle Length
w = l/2; % Vehicle width
r=0.5; % Wheel Radius
mw=m/2; % Wheel Mass

% Moments of Inertia
Ibody=m*(l^2+w^2)/12;
Jsteer=mw*r^2/4;

origin=0;

epsilon=0.05;

a=1; % Damping coefficient

K0=-1;


% Potential Field
Field = 1;
% 1 = Minima, 2 = Canyon, 3 = Sine Canyon,
% 4 = Gaussian at (xf,yf), 5 = Dual Parabaloid, 6 = Parallel Park
% 7 = Slanted Saddle (Gate), 8 = Slope

for i=2:N
    
    K=[K0];%*tanh(norm(q(1:3)-[xf;yf;thetaf]));
    
    if thetaf > pi
        thetaf = thetaf - 2*pi;
    end
    
    % Artificial Potential Function V(x,y,theta) and Derivative
    if Field == 1
        slope=1;
        V = -(1/2)*(q(1)^2 + q(2)^2); % Minima
        dV = -[slope*(q(1)-xf);slope*(q(2)-yf);0;0]; % Minima
    elseif Field == 2
        slope=1;
        V = -(q(1)^2 - slope*q(2) + q(3)^2); % Canyon
        dV = -[2*(q(1)-0);-slope*1;0;0]; % Canyon
    elseif Field == 3
        V = -((q(2)-sin(q(1)))^2 - q(1)); % Sine canyon
        dV = -[2*cos(q(1))*(q(2)-sin(q(1)))-2;2*(q(2)-sin(q(1)));0;0]; % Sine canyon
    elseif Field == 4
        sigma = norm(q0(1)-xf,q0(2)-yf);
        dV = -[(1/(2*sigma^4*pi)*(q(1)-xf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));(1/(2*sigma^4*pi)*(q(2)-yf)*exp(-((q(1)-xf)^2+(q(2)-yf)^2)/(2*sigma^2)));0;0];
    elseif Field == 5
        n=2; % Directional quadratic. Even powers.
        m=2; % Braking polynonmial
        m1 = tan(thetaf);
        m2 = tan(thetaf-pi/2);
        if thetaf==-pi/2 || pi/2
            dV = -[(n/2)*(q(1)-xf)^(n-1);(m/2)*(q(2)-yf)^(m-1);0;0];
        else
            dV = -[(n/2)*m1*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*m2*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);(n/2)*(m1*(q(1)-xf)+(q(2)-yf))^(n-1)+(m/2)*(m2*(q(1)-xf)+(q(2)-yf))^(m-1);0;0];
        end
    elseif Field == 6
        dV = -[sech(pi-1*q0(1)/2)^2;2*q(2);0;0];
    elseif Field == 7
        K=K0;
        dV = -4*[q(1)/(q(2)^2+1);(-(q(1)^2-1)*q(2)+(q(2)^4)/2+2*q(2)^2+1)/(q(2)^2+1)^2;0;0];
    elseif Field == 8
        slope=1;
        dV = -[0;slope;0;0];
    end
    
    % Dissipative Force
    
    Fd = -a*[q(5);q(6);q(7);q(8)];
    
    % Codistriburtion
    A = [-sin(q(3)+q(4))    cos(q(3)+q(4))      l*cos(q(4))     0;...
        -sin(q(3))    cos(q(3))       0            0 ];
    
    % Distribution
    D = [cos(q(4))*cos(q(3))   0 ;...
        sin(q(3))*cos(q(4))    0 ;...
        sin(q(4))/l   0 ;...
        0             1 ];
    
    G = [m 0 0 0; 0 m 0 0; 0 0 Ibody 0; 0 0 0 Jsteer];
    Gsharp = inv(G);
    
    A*D;
    
    MAX(i) = max(max(A*D));
    % Orthogonal Projections
    Pd =  (D)*inv(transpose(D)*D)*transpose(D);
    Pdperp = transpose(A)*inv(A*transpose(A))*A;
    
    %     Pd+Pdperp;
    
    % Calculate Virtual Torques
    PHI = [(dV(1)/sqrt(dV(1)^2+1));(dV(2)/sqrt(dV(2)^2+1));(dV(3)/sqrt(dV(3)^2+1));(dV(4)/sqrt(dV(4)^2+1))];
    Fstar = K*PHI;
    dellPF = [1 0 -(l)*sin(q(3))-r*sin(q(3)+q(4)) -r*sin(q(3)+q(4)) ;...
        0 1 (l)*cos(q(3))+r*cos(q(3)+q(4)) r*cos(q(3)+q(4))  ;...
        0 0 0 0 ;...
        0 0 0 0 ];...
        
    Qstar = transpose(dellPF)*Fstar;
    % Base lambda
    
    L = [((q(7)+q(8))*cos(q(3)+q(4))*(q(5))+(q(7)+q(8))*sin(q(3)+q(4))*(q(6))+l*(q(8))*sin(q(4))*(q(7)))*l^2/(1+(l^2-1)*cos(q(4))^2);...
        (q(7))*((q(5))*cos(q(3))+(q(6))*sin(q(3)))];
    
    %         Adot   = [-(q(7)+q(8))*cos(q(3)+q(4))    -(q(7)+q(8))*sin(q(3)+q(4)) -l*q(8)*sin(q(4))   0;...
    %             -(q(7))*cos(q(3))              -(q(7))*sin(q(3))           0                   0 ];
    
    %
    %     Pdperpdot=[2*l^2*cos(q(4))*(((-2*cos(q(4))*cos(q(3))^2*sin(q(3))-2*sin(q(4))*cos(q(3))^3+sin(q(4))*cos(q(3)))*(q(8))-2*(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3)+q(4))^2+2*sin(q(3)+q(4))*cos(q(3))*((cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1))*cos(q(3)+q(4))+(cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3)*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2))*cos(q(3))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), -(2*((-2*cos(q(4))*cos(q(3))^4+2*sin(q(4))*cos(q(3))^3*sin(q(3))+2*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))-2*(cos(q(3))^2-1/2)^2*cos(q(4))*(q(7)))*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/2)*sin(q(3))*cos(q(4))-sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3))*cos(q(3)+q(4))+2*cos(q(3))^2*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+2*(cos(q(3))^2-1/2)*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(4))*(q(7)))*l^2*cos(q(4))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), l*(((-4*cos(q(4))*cos(q(3))^4+4*sin(q(4))*cos(q(3))^3*sin(q(3))+3*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))+(-4*cos(q(4))*cos(q(3))^4+5*cos(q(4))*cos(q(3))^2-cos(q(4)))*(q(7)))*cos(q(3)+q(4))^3-(4*(cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/4)*sin(q(3))*cos(q(4))-(3/4)*sin(q(4))*cos(q(3)))*(q(8))+4*(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-3/4))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^2+(-(-3*cos(q(4))*cos(q(3))^3+3*sin(q(4))*cos(q(3))^2*sin(q(3))+(l^2*cos(q(4))^3+2*cos(q(4)))*cos(q(3))-l^2*cos(q(4))^2*sin(q(4))*sin(q(3)))*cos(q(3))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1)*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))+(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2))*sin(q(3)+q(4))*cos(q(3)))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
    %         -(2*((-2*cos(q(4))*cos(q(3))^4+2*sin(q(4))*cos(q(3))^3*sin(q(3))+2*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))-2*(cos(q(3))^2-1/2)^2*cos(q(4))*(q(7)))*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/2)*sin(q(3))*cos(q(4))-sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3))*cos(q(3)+q(4))+2*cos(q(3))^2*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+2*(cos(q(3))^2-1/2)*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(4))*(q(7)))*l^2*cos(q(4))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), -2*l^2*((-(2*cos(q(3))*cos(q(4))*sin(q(3))+2*cos(q(3))^2*sin(q(4))-sin(q(4)))*(cos(q(3))+1)*(cos(q(3))-1)*(q(8))-(2*cos(q(3))^2-1)*sin(q(3))*cos(q(4))*cos(q(3))*(q(7)))*cos(q(3)+q(4))^2+(2*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*(q(8))+2*(q(7))*cos(q(3))^2*cos(q(4)))*(cos(q(3))+1)*(cos(q(3))-1)*sin(q(3)+q(4))*cos(q(3)+q(4))+((cos(q(3))-1)*(cos(q(3))+1)*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2))*cos(q(3)))*cos(q(4))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), (((-4*cos(q(4))*cos(q(3))^3*sin(q(3))-4*sin(q(4))*cos(q(3))^4+3*cos(q(3))*cos(q(4))*sin(q(3))+5*cos(q(3))^2*sin(q(4))-sin(q(4)))*(q(8))+(-4*cos(q(4))*cos(q(3))^3*sin(q(3))+cos(q(3))*cos(q(4))*sin(q(3)))*(q(7)))*cos(q(3)+q(4))^3+(4*(cos(q(4))*cos(q(3))^4-sin(q(4))*cos(q(3))^3*sin(q(3))-(5/4)*cos(q(4))*cos(q(3))^2+(3/4)*sin(q(4))*cos(q(3))*sin(q(3))+(1/4)*cos(q(4)))*(q(8))+4*(q(7))*cos(q(4))*cos(q(3))^2*(cos(q(3))^2-3/4))*sin(q(3)+q(4))*cos(q(3)+q(4))^2+((3*sin(q(4))*cos(q(3))^4+3*cos(q(4))*cos(q(3))^3*sin(q(3))+(-l^2*cos(q(4))^2*sin(q(4))-3*sin(q(4)))*cos(q(3))^2-cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+2)*cos(q(3))+l^2*cos(q(4))^2*sin(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))-(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(l^2*cos(q(4))^2+cos(q(3))^2))*sin(q(3)+q(4)))*l/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
    %         l*(((-4*cos(q(4))*cos(q(3))^4+4*sin(q(4))*cos(q(3))^3*sin(q(3))+3*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))+(-4*cos(q(4))*cos(q(3))^4+5*cos(q(4))*cos(q(3))^2-cos(q(4)))*(q(7)))*cos(q(3)+q(4))^3-(4*(cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/4)*sin(q(3))*cos(q(4))-(3/4)*sin(q(4))*cos(q(3)))*(q(8))+4*(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-3/4))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^2+(-(-3*cos(q(4))*cos(q(3))^3+3*sin(q(4))*cos(q(3))^2*sin(q(3))+(l^2*cos(q(4))^3+2*cos(q(4)))*cos(q(3))-l^2*cos(q(4))^2*sin(q(4))*sin(q(3)))*cos(q(3))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1)*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))+(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2))*sin(q(3)+q(4))*cos(q(3)))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), (((-4*cos(q(4))*cos(q(3))^3*sin(q(3))-4*sin(q(4))*cos(q(3))^4+3*cos(q(3))*cos(q(4))*sin(q(3))+5*cos(q(3))^2*sin(q(4))-sin(q(4)))*(q(8))+(-4*cos(q(4))*cos(q(3))^3*sin(q(3))+cos(q(3))*cos(q(4))*sin(q(3)))*(q(7)))*cos(q(3)+q(4))^3+(4*(cos(q(4))*cos(q(3))^4-sin(q(4))*cos(q(3))^3*sin(q(3))-(5/4)*cos(q(4))*cos(q(3))^2+(3/4)*sin(q(4))*cos(q(3))*sin(q(3))+(1/4)*cos(q(4)))*(q(8))+4*(q(7))*cos(q(4))*cos(q(3))^2*(cos(q(3))^2-3/4))*sin(q(3)+q(4))*cos(q(3)+q(4))^2+((3*sin(q(4))*cos(q(3))^4+3*cos(q(4))*cos(q(3))^3*sin(q(3))+(-l^2*cos(q(4))^2*sin(q(4))-3*sin(q(4)))*cos(q(3))^2-cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+2)*cos(q(3))+l^2*cos(q(4))^2*sin(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))-(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(l^2*cos(q(4))^2+cos(q(3))^2))*sin(q(3)+q(4)))*l/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 4*(q(8))*((cos(q(3))*cos(q(4))*sin(q(3))+cos(q(3))^2*sin(q(4))-(1/2)*sin(q(4)))*cos(q(3)+q(4))^2-(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*sin(q(3)+q(4))*cos(q(3)+q(4))-(1/2)*cos(q(3))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3))))*l^2*cos(q(4))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+(8*cos(q(3))^2-4)*sin(q(3))*sin(q(3)+q(4))*cos(q(3))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
    %         0, 0, 0, 0];
    
    Pdperpdot=[2*cos(q(3))*(((-2*cos(q(4))*cos(q(3))^2*sin(q(3))-2*sin(q(4))*cos(q(3))^3+sin(q(4))*cos(q(3)))*(q(8))-2*(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3)+q(4))^2+2*cos(q(3))*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1))*cos(q(3)+q(4))+(cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3)*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2))*cos(q(4))*l^2/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), -(2*((-2*cos(q(4))*cos(q(3))^4+2*sin(q(4))*cos(q(3))^3*sin(q(3))+2*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))-2*cos(q(4))*(q(7))*(cos(q(3))^2-1/2)^2)*cos(q(3)+q(4))^2-4*cos(q(3))*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/2)*sin(q(3))*cos(q(4))-sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3)+q(4))+2*cos(q(3))^2*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+2*cos(q(4))*(q(7))*(cos(q(3))^2-1/2)*(l^2*cos(q(4))^2+cos(q(3))^2))*cos(q(4))*l^2/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), (((-4*cos(q(4))*cos(q(3))^4+4*sin(q(4))*cos(q(3))^3*sin(q(3))+3*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))+(-4*cos(q(4))*cos(q(3))^4+5*cos(q(4))*cos(q(3))^2-cos(q(4)))*(q(7)))*cos(q(3)+q(4))^3-4*cos(q(3))*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/4)*sin(q(3))*cos(q(4))-(3/4)*sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-3/4))*cos(q(3)+q(4))^2+(-cos(q(3))*(-3*cos(q(4))*cos(q(3))^3+3*sin(q(4))*cos(q(3))^2*sin(q(3))+(l^2*cos(q(4))^3+2*cos(q(4)))*cos(q(3))-l^2*cos(q(4))^2*sin(q(4))*sin(q(3)))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1)*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))+cos(q(3))*sin(q(3)+q(4))*(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)))*l/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
        -(2*((-2*cos(q(4))*cos(q(3))^4+2*sin(q(4))*cos(q(3))^3*sin(q(3))+2*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))-2*cos(q(4))*(q(7))*(cos(q(3))^2-1/2)^2)*cos(q(3)+q(4))^2-4*cos(q(3))*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/2)*sin(q(3))*cos(q(4))-sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-1/2))*cos(q(3)+q(4))+2*cos(q(3))^2*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+2*cos(q(4))*(q(7))*(cos(q(3))^2-1/2)*(l^2*cos(q(4))^2+cos(q(3))^2))*cos(q(4))*l^2/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), -(2*(-(2*cos(q(3))+2)*(cos(q(3))*cos(q(4))*sin(q(3))+cos(q(3))^2*sin(q(4))-(1/2)*sin(q(4)))*(cos(q(3))-1)*(q(8))-2*cos(q(3))*cos(q(4))*(q(7))*(cos(q(3))^2-1/2)*sin(q(3)))*cos(q(3)+q(4))^2+2*(2*cos(q(3))+2)*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*(q(8))+cos(q(3))^2*cos(q(4))*(q(7)))*(cos(q(3))-1)*cos(q(3)+q(4))+2*cos(q(3))*((cos(q(3))-1)*(cos(q(3))+1)*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)))*cos(q(4))*l^2/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), l*(((-4*cos(q(4))*cos(q(3))^3*sin(q(3))-4*sin(q(4))*cos(q(3))^4+3*cos(q(3))*cos(q(4))*sin(q(3))+5*cos(q(3))^2*sin(q(4))-sin(q(4)))*(q(8))+(-4*cos(q(4))*cos(q(3))^3*sin(q(3))+cos(q(3))*cos(q(4))*sin(q(3)))*(q(7)))*cos(q(3)+q(4))^3+4*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^4-sin(q(4))*cos(q(3))^3*sin(q(3))-(5/4)*cos(q(4))*cos(q(3))^2+(3/4)*sin(q(4))*cos(q(3))*sin(q(3))+(1/4)*cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(cos(q(3))^2-3/4))*cos(q(3)+q(4))^2+((3*sin(q(4))*cos(q(3))^4+3*cos(q(4))*cos(q(3))^3*sin(q(3))+(-l^2*cos(q(4))^2*sin(q(4))-3*sin(q(4)))*cos(q(3))^2-cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+2)*cos(q(3))+l^2*cos(q(4))^2*sin(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))-sin(q(3)+q(4))*(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(l^2*cos(q(4))^2+cos(q(3))^2)))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
        (((-4*cos(q(4))*cos(q(3))^4+4*sin(q(4))*cos(q(3))^3*sin(q(3))+3*cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3)))*(q(8))+(-4*cos(q(4))*cos(q(3))^4+5*cos(q(4))*cos(q(3))^2-cos(q(4)))*(q(7)))*cos(q(3)+q(4))^3-4*cos(q(3))*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^2*sin(q(3))+sin(q(4))*cos(q(3))^3-(1/4)*sin(q(3))*cos(q(4))-(3/4)*sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(cos(q(3))^2-3/4))*cos(q(3)+q(4))^2+(-cos(q(3))*(-3*cos(q(4))*cos(q(3))^3+3*sin(q(4))*cos(q(3))^2*sin(q(3))+(l^2*cos(q(4))^3+2*cos(q(4)))*cos(q(3))-l^2*cos(q(4))^2*sin(q(4))*sin(q(3)))*(q(8))+(q(7))*cos(q(4))*(cos(q(3))-1)*(cos(q(3))+1)*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))+cos(q(3))*sin(q(3)+q(4))*(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3)))*(q(8))+(q(7))*cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)))*l/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), l*(((-4*cos(q(4))*cos(q(3))^3*sin(q(3))-4*sin(q(4))*cos(q(3))^4+3*cos(q(3))*cos(q(4))*sin(q(3))+5*cos(q(3))^2*sin(q(4))-sin(q(4)))*(q(8))+(-4*cos(q(4))*cos(q(3))^3*sin(q(3))+cos(q(3))*cos(q(4))*sin(q(3)))*(q(7)))*cos(q(3)+q(4))^3+4*sin(q(3)+q(4))*((cos(q(4))*cos(q(3))^4-sin(q(4))*cos(q(3))^3*sin(q(3))-(5/4)*cos(q(4))*cos(q(3))^2+(3/4)*sin(q(4))*cos(q(3))*sin(q(3))+(1/4)*cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(cos(q(3))^2-3/4))*cos(q(3)+q(4))^2+((3*sin(q(4))*cos(q(3))^4+3*cos(q(4))*cos(q(3))^3*sin(q(3))+(-l^2*cos(q(4))^2*sin(q(4))-3*sin(q(4)))*cos(q(3))^2-cos(q(4))*sin(q(3))*(l^2*cos(q(4))^2+2)*cos(q(3))+l^2*cos(q(4))^2*sin(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+3*cos(q(3))^2))*cos(q(3)+q(4))-sin(q(3)+q(4))*(-(l*cos(q(4))-cos(q(3)))*(l*cos(q(4))+cos(q(3)))*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-cos(q(4)))*(q(8))+(q(7))*cos(q(4))*cos(q(3))^2*(l^2*cos(q(4))^2+cos(q(3))^2)))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), (4*(cos(q(3))*cos(q(4))*sin(q(3))+cos(q(3))^2*sin(q(4))-(1/2)*sin(q(4)))*cos(q(3)+q(4))^2-4*(cos(q(4))*cos(q(3))^2-sin(q(4))*cos(q(3))*sin(q(3))-(1/2)*cos(q(4)))*sin(q(3)+q(4))*cos(q(3)+q(4))-2*cos(q(3))*(sin(q(3))*cos(q(4))+sin(q(4))*cos(q(3))))*cos(q(4))*l^2*(q(8))/((8*cos(q(3))^4-8*cos(q(3))^2+1)*cos(q(3)+q(4))^4+8*cos(q(3))*(cos(q(3))^2-1/2)*sin(q(3))*sin(q(3)+q(4))*cos(q(3)+q(4))^3+(-8*cos(q(3))^4+(-4*l^2*cos(q(4))^2+6)*cos(q(3))^2+2*l^2*cos(q(4))^2)*cos(q(3)+q(4))^2-4*sin(q(3)+q(4))*cos(q(3))*sin(q(3))*(l^2*cos(q(4))^2+cos(q(3))^2)*cos(q(3)+q(4))+(l^2*cos(q(4))^2+cos(q(3))^2)^2), 0;...
        0, 0, 0, 0];
    
    Lambda = -Pdperpdot*q(5:8);
    
    Lambda1(:,i) = -Pdperpdot*q(5:8); % From dynamic equations
    Lambda2(:,i) = transpose(A)*L;
    
    
    Lambda = Lambda1(:,i);
    
    qdot = [q(5);q(6);q(7);q(8);Gsharp*Pd*((dV+Fd)+[0;0;0;Qstar(4)])+Gsharp*Lambda];
    
    q = q+dt*qdot;
    
    constraint(:,i) =  A*qdot(1:4); % check if non-holonomic constraint is satisfied i.e., A(q)*qdot == 0;
    
    Q(:,i) = q;
    U(:,i) = [Qstar(1)*cos(q(3))+Qstar(1)*sin(q(3));Qstar(4)];
    %     L(:,i)=Lambda;
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
    Z = 1/2*(slope*(X-xf).^2+slope*(Y-yf).^2); % Minima
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
elseif Field==8
    Z = slope*Y;
end

Zmax = max(max(abs(Z)));
Z=1*Z./Zmax;

figure(1)
set(gcf,'color','w');
clf
for i=N:1000:N
    %     clf
    subplot(2,2,[2])
    plot(t(2:i),Q(1:8,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$q$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$','$\phi$'},'Interpreter','latex','Location','best')
    %     legend({'$x$','$y$','$\theta$','$\phi$','$\dot{x}$','$\dot{y}$','$\dot{\theta}$','$\dot{\phi}$'},'Interpreter','latex','Location','best')
    grid on
    box on
    subplot(2,2,[4])
    plot(t(2:i),constraint(:,2:i),t(2:i),U(:,2:i));
    legend({'$\alpha_1 \dot{q}$','$\alpha_2 \dot{q}$','$\tau_1$','$\tau_2$'},'Interpreter','latex','Location','best'); %
    grid on
    box on
    
    subplot(2,2,[1 3])
    surf(X,Y,Z-1);
    shading interp;
    colormap jet;
    colorbar('Ticks',[min(min(Z))-1,0],...
        'TickLabels',{'Minimum','Maximum'})  ;
    %     az = 0; el = 90; view(az, el);
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    hold on
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Q(1,2:i),Q(2,2:i),'k--','Linewidth',1);
    hold on
    plotCar(Q(1,i)+l/2*cos(Q(3,i)),Q(2,i)+l/2*sin(Q(3,i)),Q(3,i),Q(4,i),l,w,r);
    
    %     viscircles([0 0],epsilon);
    %     viscircles([0 0],norm(q0(1:2)),'Color','b');
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    set(gca,'ZTick',[])
    %     axis equal
    grid off
    box on
    pause(10^-9)
    hold off
end
%%
figure(2)
plot(t,Lambda1(1,:),t,Lambda2(1,:))
