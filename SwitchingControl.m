close all;
clear Q
% Example for Lyapunov, Filipov in a Switching Scenario

% Plot the vector field (xdot in the filipov)
figure(1)
hold on
lambda=-2; 
A1 = [0 1;-1 0];
A2 = [-1 0; 0 lambda];

for x = -10:1:10
    for y=-10:1:10 
        if y>=x
            dot = A1*[x;y];
        elseif y<x
            dot = A2*[x;y];
        end
        figure(1)
        quiver(x,y,dot(1),dot(2),0.1)
        hold on
    end
end

T = 10; % simulation time
dt=0.01; % step time
N = T/dt;
t=linspace(0,T,N);

q0 = 10*randn(2,1); % Initial Conditions (make this arbitrary)

q=q0;

Q(:,1)=q0;

for i=2:N
    if q(2)>=q(1)
        dot = A1*[q(1);q(2)];
    elseif q(2)<q(1)
        dot = A2*[q(1);q(2)];
    end
    
    q=q+dt*dot;
    
    Q(:,i)=q;
end


plot(Q(1,:),Q(2,:),'--r')



