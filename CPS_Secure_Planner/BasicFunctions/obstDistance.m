function dist =  obstDistance(x, sqrtQ, v)
global xO rSafe rReact polygons
persistent verticies

if isempty(verticies)
     verticies = polygons;
end

numObst = length(verticies);
%Begin testing

 %Should do this with KDtree
dist = Inf;
if(isempty(sqrtQ))
     for i=1:numObst
     transVert = verticies{i}';
     distNow = norm(p_poly_dist(x(1),x(2), transVert(1,:), transVert(2,:)));
     %distNow = norm(sqrtQ*(-x(1:2)+xO(i,:)'))-1;
    if( distNow < dist)
     dist = distNow;
    end
     end

    
else 
    transX = sqrtQ*x(1:2);
     for i=1:numObst
         transVert = sqrtQ* verticies{i}';
         distNow = norm(p_poly_dist(transX(1),transX(2), transVert(1,:), transVert(2,:))) -1;
         %distNow = norm(sqrtQ*(-x(1:2)+xO(i,:)'))-1;
        if( distNow < dist)
         dist = distNow;
        end
     end

     dist = dist/max(abs(sqrtQ(1,1)),abs(sqrtQ(2,2)));
 
end
 

end