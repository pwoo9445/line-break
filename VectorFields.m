% Simple Vector Field for c1=1. The "simplest" case.

clear all;

T = 60; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

q0 = [7*randn(2,1);pi*randn]; % Arbitrary start position

q=q0;

Q(:,1)=q0;

ku=0.4;
kw=4;

for i=2:N
    
    POS = [q(1);q(2)];
    
    u = ku*norm(POS);
    
    phi = atan2(q(2),q(1))+pi;
    
    
    R = sqrt(q(1)^2+q(2)^2);
    F.x = -q(1)/R;
    F.y = -q(2)/R;
    dF.x.x = -q(2)^2/R^3;
    dF.x.y = q(1)*q(2)/R^3;
    dF.y.y = q(1)*q(2)/R^3;
    dF.y.x =  q(1)^2/R^3;
    
    phidot = ((dF.y.y*sin(q(3))+dF.y.x*cos(q(3)))*F.x-(dF.x.y*sin(q(3))+dF.x.x*cos(q(3)))*F.y)*u;
    
    omega = -kw*(q(3)-phi) + phidot;
    
    qdot = [cos(q(3)) 0 ; sin(q(3)) 0; 0 1]*[u;omega];
    
    q = q+dt*qdot;
    Q(:,i) = q;
    
end

j=1;
clear F dF

for x = -20:0.5:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
    for y= -20:0.5:20 
        
        phi(j) = atan2(y,x) + pi ;
        
        Fx(j) = cos(phi(j));
        Fy(j) = sin(phi(j));
        
        A = sqrt((y/x)^2+1);
        R = sqrt(x^2+y^2)^2;
        
        F.x(j) = -x/R;
        F.y(j) = -y/R;
        dF.x.x(j) = -y^2/R^3;
        dF.x.y(j) = x*y/R^3;
        dF.y.y(j) = x*y/R^3;
        dF.y.x(j) =  x^2/R^3;
        
        X(:,j) = [x;y];
        j=j+1;
    end
end

close all
figure(1)
set(gcf,'color','w');

plot(t,Q);
ylim([min(min(Q)),max(max(Q))]);
xlim([min(t),max(t)]);
legend({'$x~[m]$','$y~[m]$','$\theta~[rad]$'},'Interpreter','latex','Location','best')
box on

figure(2)
quiver(X(1,:),X(2,:),F.x,F.y,0.3)
hold on
plot(Q(1,:),Q(2,:))


%% Attractive
close all;
A0 = [10;10]

j=1;
for x = -20:0.5:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
    for y= -20:0.5:20 
       x=x-A0(1);
              y=y-A0(2);

        A = sqrt((y/x)^2+1);
        R = sqrt(x^2+y^2)^2;
        
        F.x(j) = -x/R;
        F.y(j) = -y/R;
        dF.x.x(j) = -y^2/R^3;
        dF.x.y(j) = x*y/R^3;
        dF.y.y(j) = x*y/R^3;
        dF.y.x(j) =  x^2/R^3;
        
              x=x+A0(1);
              y=y+A0(2);
        X(:,j) = [x;y];
        j=j+1;
    end
end

quiver(X(1,:),X(2,:),F.x,F.y,0.5)

%% Repulsive

clear all
close all
x0=0;
y0=0;

j=1;
for x = -20:0.5:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
    for y= -20:0.5:20 
        
        x=x-x0;
        y=y-y0;
                
        A = sqrt((y/x)^2+1);
        R = sqrt(x^2+y^2)^2;
        
        F.x(j) = x/R;
        F.y(j) = y/R;
        dF.x.x(j) = y^2/R^3;
        dF.x.y(j) = x*y/R^3;
        dF.y.y(j) = x*y/R^3;
        dF.y.x(j) = x^2/R^3;
        
        x=x+x0;
        y=y+y0;
        
           X(:,j) = [x;y];

        j=j+1;
    end
end


quiver(X(1,:),X(2,:),F.x,F.y)



%% Conservative (pass on right)
clear all;
clf
x0=10;
y0=0;

j=1;
for x = -20:0.5:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
    for y= -20:0.5:20 
        
        x=x-x0;
        y=y-y0;
                
        A = sqrt((y/x)^2+1);
        R = sqrt(x^2+y^2)^2;
        
        F.x(j) = -y/R;
        F.y(j) = x/R;
        dF.x.x(j) = y^2/R^3;
        dF.x.y(j) = x*y/R^3;
        dF.y.y(j) = x*y/R^3;
        dF.y.x(j) = x^2/R^3;
        
        x=x+x0;
        y=y+y0;
        
           X(:,j) = [x;y];

        j=j+1;
    end
end


quiver(X(1,:),X(2,:),F.x,F.y,0.6)




