%%%%Written by Alexander I. Ivanov - 2017%%%%
function dist =  obstDistance(x, sqrtQ)
%This function simply does a distance calculation via obstacle checks.
%In more complex environments this code should be improved so that each
%obstacle need not be checked. Note that this is not a true distance 
%function since when x is inside and obstacle the distance is given as 
%negative. This is necessary to help the NLP solve. 
global xO rSafe rReact polygons
persistent verticies safetyMargin

if isempty(verticies)
     verticies = polygons;
     safetyMargin = .025;
end

numObst = length(verticies);

 %Should do this with KDtree
dist = Inf;
%If we're calculating pure distance from obstacles
if(isempty(sqrtQ))
     for i=1:numObst
     transVert = verticies{i};
     inside = inpolygon(x(1), x(2), transVert(1,:), transVert(2,:));
         if(inside)
            distNow = -norm(p_poly_dist(x(1),x(2), transVert(1,:), transVert(2,:)));
         else
            distNow = norm(p_poly_dist(x(1),x(2), transVert(1,:), transVert(2,:)));
         end
     
        if( distNow < dist)
         dist = distNow;
        end
     end
     
     dist = dist-safetyMargin;
     return 

%This distance is the distance of an elipse to polygonal obstacles.
%We transform the obstacles then return the distnace to the unit circle. 
else 
     for i=1:numObst
         transVert = sqrtQ* (verticies{i} - x(1:2));
         distNow = norm(p_poly_dist(0,0, transVert(1,:), transVert(2,:))) -1;
        if( distNow < dist)
         dist = distNow;
        end
     end

     dist = dist/max(abs(sqrtQ(1,1)),abs(sqrtQ(2,2)));
 
end
 

end