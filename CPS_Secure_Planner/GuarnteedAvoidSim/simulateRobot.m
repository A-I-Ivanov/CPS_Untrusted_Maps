function [ realizedTraj, setPoint ] = simulateRobot( plannedTraj)
%This function is a simple simulator for the obstacle-detection
%and avoidance scenario. The optimal path is followed. When an unknown
%obstacle is located, the reactive controller is activated.
%Note that this simulator does not directly evaluate the safety of 
%The remaining path, but assumes that the new obstacle is obstructing
%the planned trajectory. 
global nx nu polygons rSensor thetaSensor unknownObst K velBounds
plannedControls(1,:) = plannedTraj(2+nx:nx+nu:end);
plannedControls(2,:) = plannedTraj(3+nx:nx+nu:end);
xNow = plannedTraj(2:nx+1);
realizedTraj = zeros(1);
realizedTraj(2:nx+1) = xNow;
obstacleDetected = 0;
setPoint = [];
dT = plannedTraj(1);
baseSteps =10;
reactiveSteps =baseSteps;
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
           uNow = reactiveControl(xNow,lastObstPoint, setPoint, .1,reactiveSteps, 1)
       else
           lastObstPoint = obstPoint;
            if(reactiveSteps == baseSteps)
                uNow = reactiveControl(xNow,obstPoint, setPoint, .1,reactiveSteps,0)
            else
                uNow = reactiveControl(xNow,obstPoint, setPoint, .1,reactiveSteps, 1)
            end
       end
       
       reactiveSteps = reactiveSteps-1;
   else
       uNow = plannedControls(:,i);
       dT = plannedTraj(1);
   end
   
  

   
   xNow = diffDriveKinematics(xNow, uNow, dT);
   
   xNow = checkBounds(xNow);
   
   realizedTraj(i*(nx+nu)+2:(i+1)*(nx+nu)+1) = [xNow;uNow];
   if( isempty(unknownObst) ==0 && inpolygon(xNow(1), xNow(2), unknownObst{1}(1,:),unknownObst{1}(2,:)))
       break
   end
   
    
end



    function correctedX = checkBounds(x)
        correctedX = x;
        if(x(4)<velBounds(1))
            correctedX(4) = velBounds(1);
        end
        if(x(4)>velBounds(2))
            correctedX(4) = velBounds(2);
        end
 
    end



end

