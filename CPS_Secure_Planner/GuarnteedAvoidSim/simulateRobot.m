function [ realizedTraj ] = simulateRobot( plannedTraj )
%SIMULATE Summary of this function goes here
%   Detailed explanation goes here
global nx nu polygons rSensor thetaSensor unknownObst K
plannedControls(1,:) = plannedTraj(2+nx:nx+nu:end);
plannedControls(2,:) = plannedTraj(3+nx:nx+nu:end);
xNow = plannedTraj(2:nx+1);
realizedTraj = zeros(1);
realizedTraj(2:nx+1) = xNow;
obstacleDetected = 0;
setPoint = [];
dT = plannedTraj(1);
reactiveSteps =10;
%Main simulation loop
for i =1:K-1
   obstPoint = checkObstacles(xNow, polygons, rSensor, thetaSensor);
   if(obstacleDetected || isempty(obstPoint)==0)
       if(isempty(setPoint))
            obstacleDetected =1;
            setPoint = xNow;
            lastObstPoint = obstPoint;
       end
       
       if(reactiveSteps<=1)
           break
       end
       
       if(isempty(obstPoint)) %If we can't see anything, we assume the object point remains the same
           uNow = reactiveControl(xNow,lastObstPoint, setPoint, plannedTraj(1),reactiveSteps);
       else
           lastObstPoint = obstPoint;
           uNow = reactiveControl(xNow,obstPoint, setPoint, plannedTraj(1),reactiveSteps);
       end
       
       reactiveSteps = reactiveSteps-1;
   else
       uNow = plannedControls(:,i);
       dT = plannedTraj(1);
   end
   
   xNow = diffDriveKinematics(xNow, uNow, dT);
   realizedTraj(i*(nx+nu)+2:(i+1)*(nx+nu)+1) = [xNow;uNow];
   if( inpolygon(xNow(1), xNow(2), unknownObst{1}(1,:),unknownObst{1}(2,:)))
       
       break
   end
   
    
end





end

