%%%%Written by Alexander I. Ivanov - 2017%%%%
function dist =  Dist2OccTri(x, triangle, sqrtQ)
%This function simply does a distance calculation via obstacle checks.
%In more complex environments this code should be improved so that each
%obstacle need not be checked. Note that this is not a true distance 
%function since when x is inside and obstacle the distance is given as 
%negative. This is necessary to help the NLP solve. 
global  polygons
persistent verticies

if isempty(verticies)
     verticies = polygons;
end

numObst = length(verticies);

 %Should do this with KDtree
dist = Inf;
distInside = 0;

P1.x =  triangle(1,:);
P1.y =  triangle(2,:);
%This distance is the distance of an elipse to polygonal obstacles.
%We transform the obstacles then return the distnace to the unit circle. 
 for i=1:numObst
     transVert = sqrtQ* (verticies{i} - x(1:2));
     
     %%Check for intersecting obstacles
     for k =1:length(transVert)
         if(inpolygon(transVert(1,k),  transVert(2,k), triangle(1,:), triangle(2,:)))
            tempDist =  p_poly_dist(transVert(1,k),  transVert(2,k), transVert(1,:), transVert(2,:));
            if(tempDist>distInside)
                distInside =  tempDist; %Return the maximum distance
            end
         end
     end
     
     %If we have an intersection, there is no reason to compute the
     %exterior distance
     if(distInside ==0)
         P2.x =  transVert(1,:);
         P2.y =  transVert(2,:);
         distNow = min_dist_between_two_polygons(P1,P2,false);
        if( distNow < dist)
         dist = distNow;
        end
     end
     
     
 end
 
 %%The distance inside the occlusion trinagle scaled
if(distInside >0)
   dist = -distInside/max(abs(sqrtQ(1,1)),abs(sqrtQ(2,2))); 
else
   dist = dist/max(abs(sqrtQ(1,1)),abs(sqrtQ(2,2))); %Distance outside occlusion trinagle scaled
end


end