function [X,Y] = plot_unicycle(Q,AX)
% PLOT_UNICYCLE   Function that generates lines for plotting a unicycle (bascially, an arrow).
%
%    PLOT_UNICYCLE(Q,AX) takes in the axis AX = axis(gca) and the unicycle
%    configuration Q and outputs lines (X,Y) for use, for example, as:
%       fill(X,Y,'b')
%    This must be done in the parent script file.

    x = Q(1);
    y = Q(2);
    theta = Q(3);
    
    l = 0.02*max([AX(2)-AX(1),AX(4)-AX(3)]);
    X = [x,x+l*cos(theta-2*pi/3),x+l*cos(theta),x+l*cos(theta+2*pi/3),x];
    Y = [y,y+l*sin(theta-2*pi/3),y+l*sin(theta),y+l*sin(theta+2*pi/3),y];

end