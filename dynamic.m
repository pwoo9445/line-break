%% Vector Field for Panagou Paper
clear all
T = 20; % simulation time
dt= 0.01; % step time
N = T/dt;
t=linspace(0,T,N);

q0 = [0;0;0]; % Arbitrary start position
v0 = [0;0]
Q(:,1)=q0;
V(:,1)=v0;

v=v0;
q=q0;
% System Properties

m=1;
r=1;
J=m*r^2/2;

F = cos(.2*t);
tau = sin(0.8*t);

for i=2:N
    
    
    vdot = [1/m 0 ; ...
        0 1/J]*[F(i);tau(i)];
    
    v = v + vdot*dt;
    
    qdot = [cos(q(3)) 0 ; sin(q(3)) 0; 0 1]*v;
    
    q = q+dt*qdot;
    
    
    Q(:,i) = q;
    V(:,i) = v;
    
    
end

figure(1)
clf
set(gcf,'color','w');
subplot(2,2,[1 3])
plot(Q(1,:),Q(2,:));
xlabel({'$x~[m]$'},'Interpreter','latex')
ylabel({'$y~[m]$'},'Interpreter','latex')
legend({'$\gamma$'},'Interpreter','latex','Location','best')

axis equal
subplot(2,2,[2  4])
plot(t,V(1,:),t,V(2,:),t,F,t,tau);
hold on
legend({'$v$','$\omega$','$F$','$\tau$'},'Interpreter','latex','Location','best')





