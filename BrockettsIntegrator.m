% Brockett's Integrator
clear
clf
T = 10; % simulation time
dt= 0.001; % step time
N = T/dt;
t=linspace(0,T,N);

% Prescribed Inputs
% 
% U1 = sin(0.2*t)+exp(-.02*t);
% U2 = cos(0.1*t)+exp(-0.05*t);

% Control Gains

k1=2;
k2=2;

q0=3*randn(3,1);
q=q0;
Q(:,1)=q0;

for i=2:N
    
    u1 = -k1*sign(q(1)+q(2)*q(3))*tanh(norm(q));
    u2 = -k2*sign(q(2)-q(1)*q(3))*tanh(norm(q));
    
    qdot = [1 0 ; 0 1; q(2) -q(1)]*[u1;u2];
    q = q+dt*qdot;
    
    A = [-q(2) q(1) 1];
    
    H(i) = A*qdot;
    
    Q(:,i) = q;
    U1(i)=u1;
    U2(i)=u2;
end
%%
figure(1)
set(gcf,'color','w');

for i=1:5*dt*N:N
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
    pause(0.00001)
end
%%
plot(t,H)
xlim([0 T])
ylim([-1 1])