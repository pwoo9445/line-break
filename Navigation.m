clf

x=-5:1:5;
y=-5:1:5;

x0=5;
y0=5;
theta0=pi/10;

[X,Y] = meshgrid(x,y);

Z=X.^2+Y.^2;

plot3(X,Y,Z);

OLD = [X.'; Y.'; Z.'];

NEW = [cos(theta0) -sin(theta0) 0 ;...
    sin(theta0) cos(theta0) 0;...
    0 0 1]*OLD;

hold on
plot3(NEW(1,:),NEW(2,:),NEW(3,:))

axis equal


%% Parabolic Cylinders
clf
x0=0;
y0=0;
theta0=0;
w = [0 0 ; 1 1 ; 1 2 ; 1 3 ; 2 3 ; 3 4 ];
N = 2; % Number of waypoints
w = 5*randn(N,2);
r = 2;
W = length(w); % count number of way points for nav.
track=1; % waypoint tracking index
q=[x0;y0;theta0];

for i=2:W
    
    xmin=min(w(i-1,1),w(i,1));
    xmax=max(w(i-1,1),w(i,1));
    
    dy = (w(i,2)-w(i-1,2));
    dx = (w(i,1)-w(i-1,1))
    
    if dx == 0
        x = [w(i,1) w(i,1)];
        y = [w(i-1,2) w(i,2)];
    else
        m = dy/dx;
        x = xmin:0.05:xmax;
        y = m*(x-w(i,1))+w(i,2);
        
    end
    
    [X,Y] = meshgrid(x,y);
    
    Z = (Y-(m*(X-w(i,1))+w(i,2))).^2; % Canyon
    
    Z = Z./max(max(abs(Z)));
    
    plot(x,y,w(:,1),w(:,2),'.');
    hold on
    surface(X,Y,Z)
    
end

% plot(w(:,1),w(:,2))

%% Planes

clf
w = [0 0 ; 10 10];
W = length(w); % count number of way points for nav.
r=10;
s = 1;

a = 10;
b=2;

for i=2:W
    
    xmin=min(w(i-1,1),w(i,1));
    xmax=max(w(i-1,1),w(i,1));
    
    ymin=min(w(i-1,2),w(i,2));
    ymax=max(w(i-1,2),w(i,2));
    
    x=linspace(xmin-r,xmax+r,100);
    y=linspace(ymin-r,ymax+r,100);
    
    
    dy = (w(i,2)-w(i-1,2));
    dx = (w(i,1)-w(i-1,1));
    
    m = dy/dx;
    
    alpha = atan(s);
    theta = atan2(dy,dx);
    
    [X,Y] = meshgrid(x,y);
    
    
    X = X ;
    Y = Y ;
    
        Z = -0*s*(X*cos(theta)+Y*sin(theta));

    
    Z = Z+ 10*real(sqrt(1-(((X-6).*cos(theta)+(Y-3).*sin(theta))/a).^2-(((Y-3).*cos(theta)-(X-6).*sin(theta))/b).^2));
    
    
    %     if dx == 0
    %         x = [w(i,1) w(i,1)];
    %         y = [w(i-1,2) w(i,2)];
    %     else
    %         m = dy/dx;
    %         x = xmin:.1:xmax;
    %         y = m*(x-w(i,1))+w(i,2);
    %         y = linspace(y(1)-r,y(length(y))+r,length(x-1));
    %
    %     end
    
    %     [X,Y] = meshgrid(x,y);
    %
    %     Z = sign(w(i-1,1)-w(i,1))*(Y)+sign(w(i-1,2)-w(i,2))*(m*(X-w(i,1))+w(i,2)); % Canyon
    %
    %     Z = Z./max(max(abs(Z)));
    
    plot(x,y,w(:,1),w(:,2),'.');
    hold on
    surface(X,Y,Z)
    axis equal
end

% plot(w(:,1),w(:,2))


%% Put it all together


clf
x0=0;
y0=0;
theta0=0;
l = linspace(-pi,pi,20);
w = [-l.' -sin(l).'];
N = 4; % Number of waypoints
w = 5*randn(N,2);
w(1,:)=[0 0];
r = 0.1;
W = length(w); % count number of way points for nav.
track=1; % waypoint tracking index
q=[x0;y0;theta0];

s=1; % Slope of line per unit of run

for i=2:W
    
    
    xmin=min(w(i-1,1),w(i,1));
    xmax=max(w(i-1,1),w(i,1));
    
    ymin=min(w(i-1,2),w(i,2));
    ymax=max(w(i-1,2),w(i,2));
    
    x=linspace(xmin-r,xmax+r,10);
    y=linspace(ymin-r,ymax+r,10);
    
    
    dy = (w(i,2)-w(i-1,2));
    dx = (w(i,1)-w(i-1,1));
    
    m = dy/dx;
    
    alpha = atan(s);
    theta = atan2(dy,dx);
    
    [X,Y] = meshgrid(x,y);
        
        if dx == 0
            Z = (1/2)*(X-w(i-1,1)).^2-s*(X*cos(theta)+Y*sin(theta));
            dV = -[(q(1)-w(i-1,1))-s*cos(theta);-s*sin(theta)];
            x = [w(i,1) w(i,1)];
            y = [w(i-1,2) w(i,2)];
        else
            Z = (1/2)*(Y-(m*(X-w(i,1))+w(i,2))).^2-s*(X*cos(theta)+Y*sin(theta));
            dV = -[(-m)*(q(2)-(m*(q(1)-w(i,1))))-s*cos(theta);
                (q(2)-(m*(q(1)-w(i,1))+w(i,2)))-s*sin(theta) ]; % Canyon
            x = xmin-r:0.05:xmax+r;
            y = m*(x-w(i,1))+w(i,2);
        end
    
    
    Z = Z./max(max((abs(Z))))-(i-1)*1;
    
    plot(x,y,w(:,1),w(:,2),'b.');
    hold on
    surface(X,Y,Z)
    axis equal
end

plot(w(1,1),w(1,2),'kX');


%% Gaussian

clf
l=3; % vehicle width
w=1; % vehcile length
theta=pi/2;
A=1;
x0=1;
y0=2;
sigma_x=(w*sin(theta)+l*cos(theta));
sigma_y=(w*cos(theta)+l*sin(theta));

x=linspace(-10,10,100);
y=linspace(-10,10,100);

    [X,Y] = meshgrid(x,y);

    
    Z = A*exp(-((X-x0).^2/(2*sigma_x^2)+(Y-y0).^2/(2*sigma_y^2)));
    surface(X,Y,Z)








