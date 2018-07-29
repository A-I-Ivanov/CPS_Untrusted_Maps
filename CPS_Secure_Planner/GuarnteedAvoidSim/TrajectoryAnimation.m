function [ ] = plotObstacles_Traj_simpson( trajectory,polygons, unknownObst, varargin)
% This is a plotting function which produces graphs for trajectory
% optimization plans. The solid flag is used to plot a solid line
% instead of a quiver plot
% You can also visualize the reactive set ellipses
global nx nu thetaSensor rSensor
figure
hold on

quiverScale = .3;
solid =0;
estimateReactive = 0;

trajSize = size(trajectory);

%%%%Variables for plotting
uk_obst = [];
k_obst = [];
goal = [];
start = [];
traj = [];
set_point = [];
set_speed = [];
estimateReactive =0;

if(length(varargin)>0)
    solid = varargin{1};
end

if(length(varargin)>1)
    estimateReactive = varargin{2};
end

if(length(varargin)>2)
   setPoint = varargin{3};
else
   setPoint = [];
end

if(length(varargin)>3)
   plotViews = varargin{4};
else
   plotViews = [];
end
blkSz = nx+2*nu;

timeSteps = ceil(trajSize(1)/blkSz);

[xPos, yPos, thePos, v, w, u1, u2, u1mid, u2mid] = processTraj(trajectory,blkSz);


for k=1:timeSteps

if(solid)
 
    traj = plot(xPos, yPos, 'b');
 
else
    [polR,polAng] = pol2cart(thePos,v);

    %%Plot the x,y solutions

    traj = quiver(xPos, yPos,polR*quiverScale,polAng*quiverScale,'b','AutoScale','off');

   
end

 for i = 1:length(polygons)
   k_obst = plot([polygons{i}(1,:),polygons{i}(1,1)],[ polygons{i}(2,:), polygons{i}(2,1)],'r','LineWidth',2);
    %viscircles(xO(i,:),.1) %How far away did the robot need to be from the obstacle? Plot this
 end

  for i = 1:length(unknownObst)
    uk_obst = plot([unknownObst{i}(1,:),unknownObst{i}(1,1)],[ unknownObst{i}(2,:), unknownObst{i}(2,1)],'k','LineWidth',2);
    %viscircles(xO(i,:),.1) %How far away did the robot need to be from the obstacle? Plot this

 end


if(estimateReactive)
   plotEllipses(xPos, yPos, thePos, v);
end



if(isempty(setPoint) ==0)
  %  plotObstaclePoint(setPoint);
end

if(isempty(uk_obst))
 uk_obst = plot(0,0, 'k', 'LineWidth',2 );
end


axis([0 4 0 4])
axis equal

if(plotViews)
    for i = 15:trajSize
        plotSensor(xPos(i), yPos(i), thePos(i)); 
    end
end

end


end

function [x,y,the,v,w,u1,u2,u1mid,u2mid, deltaTs] = processTraj(traj, blkSz)

sz = size(traj);

x = [];
y =[];
the = [];
v = [];
w = [];
u1 =[];
u2 = [];
u1mid = [];
u2mid =[];

deltaTs = traj(1,:);


for i=1:sz(2)
   x = [x; traj(2:(blkSz):end,i)];
   y = [y; traj(3:(blkSz):end,i)];
   the = [the; traj(4:(blkSz):end,i)];
   v = [v; traj(5:(blkSz):end,i)];
   w = [w; traj(6:(blkSz):end,i)];
   u1 = [u1; traj(7:(blkSz):end,i)];
   u2 = [u2; traj(8:(blkSz):end,i)];
   u1mid = [u1mid; traj(9:(blkSz):end,i)];
   u2mid = [u2mid; traj(10:(blkSz):end,i)];
    
end



end

function plotSensor(x,y, theta)
global nx nu thetaSensor rSensor
        rotMat = rMatrix(0,0, theta);
        triangle = [0, cos(thetaSensor), cos(thetaSensor), 0;
                    0, sin(thetaSensor), -sin(thetaSensor), 0;]*rSensor;

        triangle = rotMat(1:2,1:2)*triangle;
        triangle = [triangle(1,:) + x; triangle(2,:) + y;];
        plot(triangle(1,:), triangle(2,:), 'y','LineWidth',2);
end

function plotObstaclePoint(setPoint)
t = linspace(0,2*pi,100);
x = cos(t);
y = sin(t);


[a, Q] = interpReactive(setPoint(4));

skewPoints = chol(Q)*[x', y']';

skewPoints(1,:) = skewPoints(1,:) +a(1);
skewPoints(2,:) = skewPoints(2,:) + a(2);

rotMat = rMatrix(0,0, setPoint(3));
skewPoints = rotMat(1:2,1:2)*skewPoints;

skewPoints(1,:) = skewPoints(1,:) + setPoint(1);
skewPoints(2,:) = skewPoints(2,:) + setPoint(2);

reactive_set = plot(skewPoints(1,:), skewPoints(2,:), 'g');

[u,v] = pol2cart(setPoint(3),setPoint(4));

set_speed = quiver(setPoint(1), setPoint(2),u*quiverScale,v*quiverScale,'b','AutoScale','off', 'MaxHeadSize', 5);

set_point = plot(setPoint(1), setPoint(2), 'rx', 'MarkerSize', 8, 'LineWidth',2);


plot(setPoint(1), setPoint(2), 'rx', 'MarkerSize', 8, 'LineWidth',2);
end

function plotEllipses(x, y, the, v)
    t = linspace(0,2*pi,100);
    xC = cos(t);
    yC = sin(t);
        for i=1:length(x)
            [a, Q] = interpReactive(v(i));

            skewPoints = chol(Q)*[xC', yC']';

            skewPoints(1,:) = skewPoints(1,:) +a(1);
            skewPoints(2,:) = skewPoints(2,:) + a(2);

            rotMat = rMatrix(0,0, the(i));
            skewPoints = rotMat(1:2,1:2)*skewPoints;

            skewPoints(1,:) = skewPoints(1,:) + x(i);
            skewPoints(2,:) = skewPoints(2,:) + y(i);

            reactive_set = plot(skewPoints(1,:), skewPoints(2,:), 'g');
        end
end

