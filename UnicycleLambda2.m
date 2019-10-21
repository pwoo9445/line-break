% Vector Field for Panagou Paper
clear 
T = 40; % simulation time
dt= 0.05; % step time
N = T/dt;
t=linspace(0,T,N);

lambda=2;

q0 = [randn*5;randn*5;0]; % Arbitrary start position

qF = [randn*5;randn*5;0]; % Final Configuration

q0(3) = atan2(-q0(2)+qF(2),-q0(1)+qF(1));

q=q0;

phi = qF(3);

Q(:,1)=q0;

ku=1;
kw=4;

for i=2:N
    
    POS = [q(1)-qF(1);q(2)-qF(2)];
    
    A = [-sin(q(3)) cos(q(3)) 0];
    
    R = sqrt(q(1)^2+q(2)^2);
    F.x = (POS(1)^2-POS(2)^2)*cos(phi) + 2*POS(1)*POS(2)*sin(phi);
    F.y = (-POS(1)^2+POS(2)^2)*sin(phi) + 2*POS(1)*POS(2)*cos(phi);
    
    F.x = F.x / sqrt(F.x^2 + F.y^2);
    F.y = F.y / sqrt(F.x^2 + F.y^2);
    
    dF.x.x = 2*POS(1)*cos(phi)+2*POS(2)*sin(phi);
    dF.x.y = -2*POS(2)*cos(phi)+2*POS(1)*sin(phi);
    dF.y.y = 2*POS(2)*cos(phi)-2*POS(1)*sin(phi);
    dF.y.x =  2*POS(1)*cos(phi)+2*POS(2)*sin(phi);
    
    dF.x.x = dF.x.x / sqrt(dF.x.x^2 + dF.x.y^2);
    dF.x.y = dF.x.y / sqrt(dF.x.x^2 + dF.x.y^2);
    dF.y.x = dF.y.x / sqrt(dF.y.x^2 + dF.y.y^2);
    dF.y.y = dF.y.y / sqrt(dF.y.x^2 + dF.y.y^2);
    
    u = ku*tanh(norm(POS));
    
    phid = atan2(F.y,F.x);
    phidot = ((dF.y.y*sin(q(3))+dF.y.x*cos(q(3)))*F.x-(dF.x.y*sin(q(3))+dF.x.x*cos(q(3)))*F.y)*u;
    
    omega = -kw*(q(3)-phid) + phidot;
    
    qdot = [cos(q(3)) 0 ; sin(q(3)) 0; 0 1]*[u;omega];
    
    q = q+dt*qdot;
    Q(:,i) = q;
    
    %     Field(:,i)=[F.x;F.y];
    H(i) = A*[F.x;F.y;0];
    H0(i) = A*qdot;
    PHID(i)=phid;
end

j=1;
for x = -20:.5:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
    for y= -20:.5:20
        
        POS = [x-qF(1);y-qF(2)];
        
        F.x = (POS(1)^2-POS(2)^2)*cos(phi) + 2*POS(1)*POS(2)*sin(phi);
        F.y = (-POS(1)^2+POS(2)^2)*sin(phi) + 2*POS(1)*POS(2)*cos(phi);
        
        F.x = F.x / sqrt(F.x^2 + F.y^2);
        F.y = F.y / sqrt(F.x^2 + F.y^2);
        
        Field(:,j)=[F.x;F.y];
        X(:,j) = [x;y];
        
        j=j+1;
    end
end


figure(1)
clf
set(gcf,'color','w');
subplot(2,2,[1 3])
plot(Q(1,:),Q(2,:),q0(1),q0(2),'kO',qF(1),qF(2),'kX');
hold on
quiver(X(1,:),X(2,:),Field(1,:),Field(2,:),1);
xlim([-20 20])
ylim([-20 20])
xlabel({'$x~[m]$'},'Interpreter','latex')
ylabel({'$y~[m]$'},'Interpreter','latex')

axis equal
subplot(2,2,[2])
plot(t,Q);
hold on
plot(t,ones(3,N).*qF,'k--')
ylim([min(min(Q)),max(max(Q))]);
xlim([min(t),max(t)]);
legend({'$x~[m]$','$y~[m]$','$\theta~[rad]$'},'Interpreter','latex','Location','best')
box on
subplot(2,2,[4])
plot(t,H,t,H0,t,PHID,t,Q(3,:))
legend({'$H=A(q)F(q)$','$A(q)\dot{q}$','$\phi_d$','$\theta$'},'Interpreter','latex','Location','best')





