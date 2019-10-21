function plotDisk(X,Y,THETA,PHI,RADIUS)

% Generate 3D Disk at (X,Y)

theta=-pi:pi/10:pi;
x=RADIUS*cos(theta);
y=zeros(1,length(theta));
z=RADIUS*sin(theta)+RADIUS;

% Rotate Disk by angle

Xnew = X+x*cos(THETA);
Ynew = Y+x*sin(THETA);

% Generate PHI indicator

% Horizontal Line

LineStart = [X;Y;RADIUS];

Segment = RADIUS*[1; 0; 0];
% Rotate by PHI
Segment = [cos(PHI) 0 -sin(PHI) ; 0 1 0 ; sin(PHI) 0 cos(PHI)]*Segment;
% Rotate by THETA
Segment = [cos(THETA) -sin(THETA) 0 ; sin(THETA) cos(THETA) 0 ; 0 0 1]*Segment;

LineEnd = LineStart+Segment;
Line = [LineStart LineEnd];

% Plot Disk
plot3(Xnew,Ynew,z,'k','Linewidth',2)
hold on
% Plot PHI line
line(Line(1,:),Line(2,:),Line(3,:),'LineWidth',1,'color','k')

