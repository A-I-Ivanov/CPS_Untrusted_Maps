function dist =  obstDistance(x, sqrtQ, v)
global xO rSafe rReact polygons
persistent verticies

if isempty(verticies)
     verticies = polygons;
end

numObst = size(xO);
%Begin testing
numObst = 2;
 %Should do this with KDtree
dist = Inf;
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
 
 
 %eigs(sqrtQ)
end