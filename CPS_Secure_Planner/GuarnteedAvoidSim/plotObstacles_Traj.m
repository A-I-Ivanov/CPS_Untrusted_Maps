function [ ] = plotObstacles_Traj( trajectory,polygons, unknownObst )
%PLOTOBSTACLES_TRAJ Summary of this function goes here
%   Detailed explanation goes here
global nx nu
figure
hold on
[u,v] = pol2cart(trajectory(4:(nx+nu):end),trajectory(5:(nx+nu):end));

%%Plot the x,y solutions

quiver(trajectory(2:(nx+nu):end), trajectory(3:(nx+nu):end),u,v,'b');
axis([-1 2 -1 2])



for i = 1:length(unknownObst)
plot([polygons{i}(1,:),polygons{i}(1,1)],[ polygons{i}(2,:), polygons{i}(2,1)],'r')
plot([unknownObst{i}(1,:),unknownObst{i}(1,1)],[ unknownObst{i}(2,:), unknownObst{i}(2,1)],'r')
%viscircles(xO(i,:),.1) %How far away did the robot need to be from the obstacle? Plot this

end

end

