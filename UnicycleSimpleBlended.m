%% Vector Field for simple attractive vector field
clear
T = 30; % simulation time
dt= 0.05; % step time
N = T/dt;
t=linspace(0,T,N);

NA=2;
AgentRad=1;
q0 = [10*randn(2,NA);mod(randn(1,NA)*pi,pi)]; % Arbitrary Start position (location of repulsive vector field)
qF = [10*randn(2,NA);mod(randn(1,NA)*pi,pi)]; % Arbitrary Final position (location of repulsive vector field)
RA = 2*AgentRad*ones(1,NA);
mindistance.A = zeros(NA,NA);

nR = 5; % Number of repulsive points
qR = 10*randn(2,nR); % Location of Conservative Vector Field
RR = 1*ones(1,nR);
mindistance.R = zeros(NA,nR);

nO = 20; % Number of conservative points
qO = 10*randn(2,nO); % Location of Conservative Vector Field
RO = 1*ones(1,nO); % Maximum Radius of contribution
mindistance.O = zeros(NA,nO);

Q = zeros(3,NA,N);

Q(:,:,1)=q0;

ku=1.5;
kw=5;

for i=2:N
    for kk=1:NA
        q = Q(:,kk,i-1);
        clear F dF sigma
        A = [-sin(q(3)) cos(q(3)) 0];
        % Attractive Vector Field
        POS = [q(1)-qF(1,kk);q(2)-qF(2,kk)];
        if norm(POS) < 0.1
            u=0; % Stop condition
        else
            u = ku*tanh(norm(POS)); % Control input
        end
        
        R = sqrt(POS(1)^2+POS(2)^2);
        F.x = -POS(1)/R;
        F.y = -POS(2)/R;
        dF.x.x = -POS(2)^2/R^3;
        dF.x.y = POS(1)*POS(2)/R^3;
        dF.y.y = POS(1)*POS(2)/R^3;
        dF.y.x =  POS(1)^2/R^3;
        
        %  Repulsive Vector Field
        for j=1:nR
            POS = [q(1)-qR(1,j);q(2)-qR(2,j)];
            R = sqrt(POS(1)^2+POS(2)^2);
            if R-RR(j) < 0
                sigma=1;
            else
                sigma = exp(-(R-RR(j))^2 / ((1)*RR(j)^2));
            end
            F.x = F.x*(1-sigma) + sigma*POS(1)/R;
            F.y = F.y*(1-sigma) + sigma*POS(2)/R;
            dF.x.x = dF.x.x*(1-sigma) + sigma*POS(2)^2/R^3;
            dF.x.y = dF.x.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
            dF.y.y = dF.y.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
            dF.y.x = dF.y.x*(1-sigma) + sigma*POS(1)^2/R^3;
        end
        
        % Conservative Vector Field
        for j=1:nO
            POS = [q(1)-qO(1,j);q(2)-qO(2,j)];
            R = sqrt(POS(1)^2+POS(2)^2);
            if R-RO(j) < 0
                sigma=1;
            else
                sigma = exp(-(R-RO(j))^2 / (2*RO(j)^2));
            end
            
            F.x = F.x*(1-sigma) - sigma*POS(2)/R;
            F.y = F.y*(1-sigma) + sigma*POS(1)/R;
            dF.x.x = dF.x.x*(1-sigma) + sigma*POS(2)^2/R^3;
            dF.x.y = dF.x.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
            dF.y.y = dF.y.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
            dF.y.x = dF.y.x*(1-sigma) + sigma*POS(1)^2/R^3;
        end
        
        % Conservative Vector Field for other agents
        for j=1:NA
            if j~=kk
                POS = [q(1)-Q(1,j,i-1);q(2)-Q(2,j,i-1)];
                R = sqrt(POS(1)^2+POS(2)^2);
                if R < RA(j)/2
                          if R-RA(j) < 0
                sigma=1;
            else
                    sigma = exp(-(R-RA(k))^2 / ((1/8)*RA(j)^2));
            end
                    
                    
                    
                    F.x = F.x*(1-sigma) + sigma*POS(1)/R;
                    F.y = F.y*(1-sigma) + sigma*POS(2)/R;
                    dF.x.x = dF.x.x*(1-sigma) + sigma*POS(2)^2/R^3;
                    dF.x.y = dF.x.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
                    dF.y.y = dF.y.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
                    dF.y.x = dF.y.x*(1-sigma) + sigma*POS(1)^2/R^3;
                    
                elseif R < RA(kk)
                              if R < RA(j)/2
                          if R-RA(j) < 0
                sigma=1;
            else
                    sigma = exp(-(R-RA(j))^2 / ((1/2)*RA(j)^2));
            end
                    
                    
                    
                    F.x = F.x*(1-sigma) - sigma*POS(2)/R;
                    F.y = F.y*(1-sigma) + sigma*POS(1)/R;
                    dF.x.x = dF.x.x*(1-sigma) + sigma*POS(2)^2/R^3;
                    dF.x.y = dF.x.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
                    dF.y.y = dF.y.y*(1-sigma) + sigma*POS(1)*POS(2)/R^3;
                    dF.y.x = dF.y.x*(1-sigma) + sigma*POS(1)^2/R^3;
                end
                
            end
        end
        
        F.x = F.x / sqrt(F.x^2 + F.y^2);
        F.y = F.y / sqrt(F.x^2 + F.y^2);
        dF.x.x = dF.x.x / sqrt(dF.x.x^2 + dF.x.y^2);
        dF.x.y = dF.x.y / sqrt(dF.x.x^2 + dF.x.y^2);
        dF.y.x = dF.y.x / sqrt(dF.y.x^2 + dF.y.y^2);
        dF.y.y = dF.y.y / sqrt(dF.y.x^2 + dF.y.y^2);
        
        phid = atan2(F.y,F.x);
        phidot = ((dF.y.y*sin(q(3))+dF.y.x*cos(q(3)))*F.x-(dF.x.y*sin(q(3))+dF.x.x*cos(q(3)))*F.y)*u;
        
        omega = -kw*(q(3)-phid) + phidot;
        
        u/omega
        
        if norm(POS) < 0.1
            omega=0; % Stop condition
        else
        
        qdot = [cos(q(3)) 0 ; sin(q(3)) 0; 0 1]*[u;omega];
        H(kk,i)=A*qdot;
        h(kk,i)=A*[F.x;F.y;0];
        
        
        q = q+dt*qdot;
        Q(:,kk,i) = q;
        
        %         %     Field(:,i)=[F.x;F.y];
        %         H(i) = A*[F.x;F.y;0];
        %         H0(i) = A*qdot;
        %         PHID(i)=phid;
    end
    end
    end
end


for i=1:N
    
    X(:,i)=Q(1,:,i)';
    Y(:,i)=Q(2,:,i)';
    THETA(:,i)=Q(3,:,i)';
end


% clear F
% j=1;
% for x = -20:0.2:20 %:(max(Q(1,:))-min(Q(1,:)))/20:max(Q(1,:))
%     for y= -20:0.2:20
%
%         % Attractive Vector Field
%         POS = [x-qF(1);y-qF(2)];
%         R = sqrt(POS(1)^2+POS(2)^2);
%         F.x = -POS(1)/R;
%         F.y = -POS(2)/R;
%
%         %         % Repulsive Vector Field
%         for k=1:nR
%
%             POS = [x-qR(1);y-qR(2)];
%             R = sqrt(POS(1)^2+POS(2)^2);
%             sigma = exp(-R^2 / (1*RR(k)^2));
%             F.x = F.x*(1-sigma) + sigma*POS(1)/R;
%             F.y = F.y*(1-sigma) + sigma*POS(2)/R;
%         end
%
%         % Conservative Vector Field
%         for k=1:nO
%
%             POS = [x-qO(1,k);y-qO(2,k)];
%             R = sqrt(POS(1)^2+POS(2)^2);
%             sigma = exp(-R^2 / (4*RO(k)^2));
%             F.x = F.x*(1-sigma) - sigma*POS(2)/R;
%             F.y = F.y*(1-sigma) + sigma*POS(1)/R;
%         end
%
%         F.x = F.x / sqrt(F.x^2 + F.y^2);
%         F.y = F.y / sqrt(F.x^2 + F.y^2);
%
%         Field(:,j)=[F.x;F.y];
%         X(:,j) = [x;y];
%         j=j+1;
%     end
% end
%%
figure(1)
set(gcf,'color','w');
for i=2:max(5,NA):N
    clf
    plot(q0(1,:),q0(2,:),'kO',qF(1,:),qF(2,:),'kX');
    hold on
    for j=1:nO
        viscircles(qO(:,j)',RO(j),'color','r');
    end
    for j=1:nR
        viscircles(qR(:,j)',RR(j),'color','b');
    end
    hold on
    for j =1:NA
        plot(X(j,2:i),Y(j,2:i),'--');
        hold on
        plotAgent(X(j,i),Y(j,i),THETA(j,i),AgentRad,0.5)
        viscircles([X(j,i) Y(j,i)],RA(j),'color','g');
    end
    
    %     xlim([min(min(Q(1,:,:))) max(max(Q(1,:,:)))])
    %     ylim([min(min(Q(2,:,:))) max(max(Q(2,:,:)))])
    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    
    axis equal
    
    pause(10^-6)
    
end
% subplot(2,2,[2])
% plot(t,Q);
% hold on
% plot(t,ones(3,N).*qF,'k--')
% ylim([min(min(Q)),max(max(Q))]);
% xlim([min(t),max(t)]);
% legend({'$x~[m]$','$y~[m]$','$\theta~[rad]$'},'Interpreter','latex','Location','best')
% box on
% subplot(2,2,[4])
% plot(t,H,t,H0,t,PHID,t,Q(3,:))
% legend({'$H=A(q)F(q)$','$A(q)\dot{q}$','$\phi_d$','$\theta$'},'Interpreter','latex','Location','best')





