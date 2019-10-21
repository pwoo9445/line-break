clf;
p=[1;0;0];
lambda=1;
for x=-10:2:10
    for y=-10:2:10
        for z=-10:2:10
            
        r = [x;y;z];
        F = lambda*dot(p,r)*r;
        
        F = F/norm(F);
        figure(1)
        quiver3(r(1),r(2),r(3),F(1),F(2),F(3));
        hold on
        pause(0.0000001)
        
        end
    end
end


 