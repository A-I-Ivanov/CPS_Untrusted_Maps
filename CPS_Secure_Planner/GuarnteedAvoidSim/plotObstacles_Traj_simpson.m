function [ ] = plotObstacles_Traj_simpson( trajectory,polygons, unknownObst, xStart,xT, varargin)
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
if(solid)
 traj = plot(trajectory(2:(blkSz):end), trajectory(3:(blkSz):end), 'b');

else
    [u,v] = pol2cart(trajectory(4:(blkSz):end),trajectory(5:(blkSz):end));

    %%Plot the x,y solutions

    traj = quiver(trajectory(2:(blkSz):end), trajectory(3:(blkSz):end),u*quiverScale,v*quiverScale,'b','AutoScale','off');

   
end

 for i = 1:length(polygons)
   k_obst = plot([polygons{i}(1,:),polygons{i}(1,1)],[ polygons{i}(2,:), polygons{i}(2,1)],'r','LineWidth',2);
    %viscircles(xO(i,:),.1) %How far away did the robot need to be from the obstacle? Plot this

 end

  for i = 1:length(unknownObst)
    uk_obst = plot([unknownObst{i}(1,:),unknownObst{i}(1,1)],[ unknownObst{i}(2,:), unknownObst{i}(2,1)],'k','LineWidth',2);
    %viscircles(xO(i,:),.1) %How far away did the robot need to be from the obstacle? Plot this

 end

 

sV = trajectory(5:(nx+2*nu):end);
sX = trajectory(2:(nx+2*nu):end);
sY = trajectory(3:(nx+2*nu):end);
sThe = trajectory(4:(nx+2*nu):end);
if(estimateReactive)
    t = linspace(0,2*pi,100);
    x = cos(t);
    y = sin(t);
        for i=1:length(sX)
            [a, Q] = interpReactive(sV(i));

            skewPoints = chol(Q)*[x', y']';

            skewPoints(1,:) = skewPoints(1,:) +a(1);
            skewPoints(2,:) = skewPoints(2,:) + a(2);

            rotMat = rMatrix(0,0, sThe(i));
            skewPoints = rotMat(1:2,1:2)*skewPoints;

            skewPoints(1,:) = skewPoints(1,:) + sX(i);
            skewPoints(2,:) = skewPoints(2,:) + sY(i);

            reactive_set = plot(skewPoints(1,:), skewPoints(2,:), 'g');
        end
end

if(isempty(setPoint) ==0)
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



goal = plot(xT(1), xT(2), 'og','markers',8, 'MarkerFaceColor','g');
start = plot(xStart(1), xStart(2), 'og','markers',8, 'MarkerFaceColor','b');

if(isempty(uk_obst))
 uk_obst = plot(0,0, 'k', 'LineWidth',2 );
end

if(isempty(setPoint))
    if(estimateReactive)
        legend([traj k_obst uk_obst reactive_set goal start],{'Trajectory','Known Obst.','Unknown Obst.', 'Reactive Bound', 'Goal' , 'Start'} ,'Location','Northwest','FontSize',14)
    else 
        legend([traj k_obst uk_obst goal start],{'Trajectory','Known Obst.','Unknown Obst.','Goal' , 'Start'} ,'Location','Northwest','FontSize',14)
    end
else
    legend([traj k_obst uk_obst reactive_set set_speed set_point goal start], ...
           {'Trajectory','Known Obst.','Unknown Obst.','Reactive Bound' , 'Detection Speed','Detection Point' ,'Goal', 'Start'},'Location','Northwest','FontSize',14)
end

xlabel('X in meters', 'FontSize',16)
ylabel('Y in meters', 'FontSize',16)



axis([0 4 0 4])
axis equal

if(plotViews)
    for i = 5:10
        rotMat = rMatrix(0,0, sThe(i));
        triangle = [0, cos(thetaSensor), cos(thetaSensor), 0;
                    0, sin(thetaSensor), -sin(thetaSensor), 0;]*rSensor;

        triangle = rotMat(1:2,1:2)*triangle;
        triangle = [triangle(1,:) + sX(i); triangle(2,:) + sY(i);];
        plot(triangle(1,:), triangle(2,:), 'y','LineWidth',2);
    end
end

end

