%%%%Written by Alexander I. Ivanov - 2017%%%%
function dist =  Dist2OccTri(triangle)
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

P1.x =  triangle(1,:);
P1.y =  triangle(2,:);
%This distance is the distance of an elipse to polygonal obstacles.
%We transform the obstacles then return the distnace to the unit circle. 

 for i=1:numObst
    transVert = ([verticies{i}(1,:); verticies{i}(2,:)]);
    plot([transVert(1,:),transVert(1,1)] , [transVert(2,:),transVert(2,1)], 'y');
     %%Check for intersecting obstacles
    P2.x =  [transVert(1,:), transVert(1,1)];
    P2.y =   [transVert(2,:), transVert(2,1)];
    distNow = min_dist_between_two_polygons(P1,P2,false);
  
    if(distNow ==0)
        
        [intersectionX,intersectionY]  = polybool('intersection',P1.x, P1.y,P2.x, P2.y);
        distNow = - polyarea(intersectionX,intersectionY); %use intersectional area as distance proxy

        
        if(distNow>0)
           notGood =1; 
        end
    end
    if( distNow < dist)
     dist = distNow;
    end
    
    

 end
 

end