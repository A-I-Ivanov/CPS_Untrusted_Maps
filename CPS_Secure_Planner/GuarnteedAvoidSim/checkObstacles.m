function [ obstPoint ] = checkObstacles (xNow, polygons, rSensor, thetaSensor)
global unknownObst;
%CHECKOBSTACLES: This function does ray tracing to perform a visibility
%check
%   Detailed explanation goes here

%The unknown obstacles are not in my map. We create and adversarial
%obstacle which is just around the corner from the known obstacle
persistent numRays verticies rays
if(isempty(verticies))
    numRays = floor(radtodeg(thetaSensor));
    %We generate a triangle looking east from the origin
    verticies = zeros(2,4);
    halfTheta = thetaSensor/2;
    verticies(:,2)= [cos(halfTheta)*rSensor; sin(halfTheta)*rSensor;];
    verticies(:,3)= [cos(-halfTheta)*rSensor; sin(-halfTheta)*rSensor;];
    
    rays = zeros( 2, numRays);
    degInRad = 2*pi/360;
    thetaNow = -halfTheta;
    for i=1:numRays
        rays(:,i) = [cos(thetaNow)*rSensor; sin(thetaNow)*rSensor;];
        thetaNow = thetaNow+degInRad;
    end
    
end
obstPoint = [];
x3d = [xNow(1:2);0];
minDist = inf; 
rotRays =[];
%%This is the main function. We first intersect unknown polygons with 
%%our FOV. If no intersection exists, we move on. If its in our FOV
%%we do ray tracing for known visible obstacles to determine the 
%%closest unknown obstacle point. 
for i=1:length(unknownObst)
   rotMat = rMatrix(0,0,xNow(3));
   rotVert = rotMat(1:2,1:2)*verticies;
   rotVert(1,:) = rotVert(1,:)+xNow(1);
   rotVert(2,:) = rotVert(2,:)+xNow(2);
   %plot(rotVert(1,:), rotVert(2,:),'b')
   
   %Can we see the unknown obstacle?
   [x_u,y_u, un_obst_inx] = polyxpoly(rotVert(1,:), rotVert(2,:), unknownObst{i}(1,:), unknownObst{i}(2,:));
   if(isempty(x_u) && inpolygon(unknownObst{i}(1,1), unknownObst{i}(2,1),rotVert(1,:), rotVert(2,:)) ==0)
       continue; %Do nothing if we cant see the unknown obstacle
   else 
       
       if(isempty(rotRays)) %Only calculate rotated rays if we need them 
           rotRays = rotMat(1:2,1:2)*rays; %rotate rays
           rotRays(1,:) = rotRays(1,:)+xNow(1); %translate Rays
           rotRays(2,:) = rotRays(2,:)+xNow(2);
          
       end
       
       
       visiPoly = zeros(1,length(polygons)); %The visible polygons
       for j =1:length(polygons)
           [x_k, y_k, ~] = polyxpoly(rotVert(1,:), rotVert(2,:), polygons{j}(1,:), polygons{j}(2,:));
           if(isempty(x_k))
               continue
           else
               visiPoly(j) = j; %We can see the jth polygon
           end
       end
       
       if(sum(visiPoly) ==0)
           for j =1:numRays
               [x_u,y_u, ~] = polyxpoly([xNow(1) rotRays(1,j)], [xNow(2) rotRays(2,j) ], unknownObst{i}(1,:), unknownObst{i}(2,:));
               [d2, pt] = calcMinDist(x_u, y_u, xNow(1:2));

               %If the distance to the unknown obstacle is less than that to
               %the known obstacle, the unknown obstacle is not occulded
               if(d2<minDist) 
                  obstPoint = pt;
                  minDist = d2;
               end
           end
       else 
         for k = find(visiPoly) %Ray trace for each visible known obstacle
           for j =1:numRays
               [x_k,y_k, ~] = polyxpoly([xNow(1) rotRays(1,j)], [xNow(2) rotRays(2,j) ], polygons{k}(1,:), polygons{k}(2,:));
               [x_u,y_u, ~] = polyxpoly([xNow(1) rotRays(1,j)], [xNow(2) rotRays(2,j) ], unknownObst{i}(1,:), unknownObst{i}(2,:));
               [d1, ~] = calcMinDist(x_k, y_k, xNow(1:2));
               [d2, pt] = calcMinDist(x_u, y_u, xNow(1:2));

               %If the distance to the unknown obstacle is less than that to
               %the known obstacle, the unknown obstacle is not occulded
               if(d1>d2 && d2<minDist) 
                  obstPoint = pt;
                  minDist = d2;
               end
           
           end
           
          end
       end
     
   end
       
end



%%%Begin helper functions 

    function [dist, point] = calcMinDist(x, y, pose)
        dist = inf;
        if(isempty(x))
           %we didn't see anything
           point = [inf; inf;];
           return; 
        end
        for u = 1:length(x)
            temp = norm(pose-[x(u);y(u)]);
            if(temp<dist)
                dist=temp;
                point = [x(u);y(u);];
            end
        end
    end

  
    
end



