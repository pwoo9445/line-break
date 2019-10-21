% Brockett's Integrator
clear
clf
T = 10; % simulation time
dt= 0.0001; % step time
N = T/dt;
t=linspace(0,T,N);

% Prescribed Inputs
%
% U1 = sin(0.2*t)+exp(-.02*t);
% U2 = cos(0.1*t)+exp(-0.05*t);

% Control Gains

k1=2;
k2=2;
k=1;
k3=2;

q0=3*randn(3,1);
q=q0;
Q(:,1)=q0;

% Coordinate Change
phi = pi/6;

Rhat = [cos(phi) sin(phi) ;...
    -sin(phi) cos(phi)];

R = [Rhat zeros(2,1)
    zeros(1,2) 1];

z0 = R*q0;


ztilde=z0;

Z(:,1)=z0;
Ztilde(:,1)=ztilde;

for i=2:N
    
    % Q coordinates
    if abs(q(1)) < 10^(-6)
        
        u1 = -k3*q(3);
        u2 = -k2*q(2);
    else
        u1 = -k1*q(1);
        u2 = (k*q(3)-(q(3)-q(1)*q(2))*u1/q(1))/q(1) ;
        
    end
    
    qdot = [1 0 ; 0 1; q(2) -q(1)]*[u1;u2];
    q = q+dt*qdot;
    
    F = [q(1)^2;q(1)*q(2);q(1)*q(3)];
    
    E(i) = -dot(qdot,F);
    tau(:,i) = cross(qdot,F);

    
    % Z coordinates
    
    z = R*q;
    V = Rhat*[u1;u2];
    
    zdot = [1 0 ; 0 1; z(2) -z(1)]*V;
    z = z+dt*zdot;
    
    % z tilde coordinates
    
    if abs(ztilde(1)/ztilde(2)-tan(phi)) < 10^(-6)
        
        vtilde1 = V(1);%-k3*z(3);
        vtilde2 = V(2);%-k2*(cos(phi)*z(1)+sin(phi)*z(2));
    else
        vtilde1 = V(1);%-k1*(cos(phi)*z(2)-sin(phi)*z(1));
        vtilde2 = (k*q(3)-(q(3)-q(1)*q(2))*vtilde1/q(1))/q(1) ;
        
    end
    
    Vtilde = [vtilde1;vtilde2];
    
    ztildedot = [1 0 ; 0 1; ztilde(2) -ztilde(1)]*Vtilde;
    ztilde = ztilde+dt*ztildedot;

    
    Q(:,i) = q;
    Z(:,i) = z;
    Ztilde(:,i)=ztilde;
    U1(i)=u1;
    U2(i)=u2;
    
    V1(i)=V(1);
    V2(i)=V(2);
    
    Vtilde(:,i) = Vtilde;
end
%%
figure(1)
set(gcf,'color','w');

figure(2)
set(gcf,'color','w');

for i=1:10*dt*N:N
    figure(1)
    clf
    subplot(2,2,[1 3])
    plot3(Q(1,1:i),Q(2,1:i),Q(3,1:i),'--',Q(1,i),Q(2,i),Q(3,i),'o');
    xlim([min(Q(1,:)) max(Q(1,:))])
    ylim([min(Q(2,:)) max(Q(2,:))])
    zlim([min(Q(3,:)) max(Q(3,:))])
    box on
    grid on
    subplot(2,2,2)
    plot(t(1:i),Q(:,1:i))
    xlim([0 T])
    ylim([min(min(Q)) max(max(Q))])
    legend({'$x_1$','$x_2$','$x_3$'},'Interpreter','latex','Location','best')
    box on
    subplot(2,2,4)
    plot(t(1:i),U1(1:i),t(1:i),U2(1:i))
    xlim([0 T])
    ylim([min(min(U1),min(U2)) max(max(U1),max(U2))])
    legend({'$u_1$','$u_2$'},'Interpreter','latex','Location','best')
    box on
    
    figure(2)
    clf
    subplot(2,2,[1 3])
    plot3(Z(1,1:i),Z(2,1:i),Z(3,1:i),'--',Z(1,i),Z(2,i),Z(3,i),'o');
    xlim([min(Z(1,:)) max(Z(1,:))])
    ylim([min(Z(2,:)) max(Z(2,:))])
    zlim([min(Z(3,:)) max(Z(3,:))])
    box on
    grid on
    subplot(2,2,2)
    plot(t(1:i),Z(:,1:i))
    xlim([0 T])
    ylim([min(min(Z)) max(max(Z))])
    legend({'$z_1$','$z_2$','$z_3$'},'Interpreter','latex','Location','best')
    box on
    subplot(2,2,4)
    plot(t(1:i),V1(1:i),t(1:i),V2(1:i))
    xlim([0 T])
    ylim([min(min(V1),min(V2)) max(max(V1),max(V2))])
    legend({'$v_1$','$v_2$'},'Interpreter','latex','Location','best')
    box on
    pause(0.00000001)
end


%%

plot(t,Ztilde,t,Z,t,Q)
legend({'$\tilde{z}_1$','$\tilde{z}_2$','$\tilde{z}_3$','$z_1$','$z_2$','$z_3$','$q_1$','$q_2$','$q_3$'},'Interpreter','latex','Location','best')

% plot(t,H)
% xlim([0 T])
% ylim([-1 1])