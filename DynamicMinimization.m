
clear Z
F=P.*Q(6,:)-Q(1,:).*cos(Q(3,:))-Q(2,:).*sin(Q(3,:));
tau=-P.*(Q(4,:).*cos(Q(3,:))+Q(5,:).*sin(Q(3,:)))-Q(3,:);


Z(:,1)=q0;
z = q0;

for i=2:N
    
    lambda=z(6)*(z(4)*cos(z(3))+z(5)*sin(z(3))); % lambda from non-holonomic ezuations
    %     if F==0
    %         F=0;
    %     else
    %     F = F/norm(F);
    %     end
    %     if tau==0
    %         tau=0;
    %     else
    %     tau = tau/norm(tau);
    %     end
    
    u = [F(i);tau(i)];
    
    zdot = [z(4) ; z(5) ; z(6) ; 1/m*(F(i)*cos(z(3))-lambda*sin(z(3))) ; 1/m*(F(i)*sin(z(3))+lambda*cos(z(3))) ; tau(i) / J];
    
    %         figure(2)
    %         plot(t(i),z(4)*sin(z(3))+z(5)*cos(z(3)),'.')
    %         hold on
    
    z = z+dt*zdot;
    Z(:,i) = z;
    U(:,i) = u;
    P(i)=p;
    
end

figure(3)
set(gcf,'color','w');

for i=N
    clf
    subplot(2,2,[2])
    plot(t(2:i),Z(1:3,2:i));
    xlabel({'$t~[s]$'},'Interpreter','latex')
    ylabel({'$z$'},'Interpreter','latex')
    legend({'$x$','$y$','$\theta$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[4])
    plot(t(2:i),U(1,2:i),t(2:i),U(2,2:i));
    legend({'$F$','$\tau$'},'Interpreter','latex','Location','best')
    grid on
    
    subplot(2,2,[1 3])
    plot(q0(1,:),q0(2,:),'kO');
    hold on
    plot(Z(1,2:i),Z(2,2:i),'--');
    hold on
%     plotAgent(Z(1,i),Z(2,i),Z(3,i),r,r/2)
    viscircles([0 0],epsilon)
        viscircles([0 0],1);
        viscircles([0 0],1,'Color','b');

    xlabel({'$x~[m]$'},'Interpreter','latex')
    ylabel({'$y~[m]$'},'Interpreter','latex')
    axis equal
    pause(10^-6)
    hold off
end


