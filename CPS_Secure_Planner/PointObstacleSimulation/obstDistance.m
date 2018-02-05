function dist =  obstDistance(x, sqrtQ, v)
global xO rSafe
persistent verticies

if isempty(verticies)
    
end

numObst = size(xO);
 %Should do this with KDtree
 dist = Inf;
 
 for i=1:numObst(1)
     distNow = norm(sqrtQ*(-x(1:2)+xO(i,:)'))-1;
    if( distNow < dist)
     dist = distNow;
    end
 end
 
 dist = dist/max(abs(sqrtQ(1,1)),abs(sqrtQ(2,2)));
 
 
 %eigs(sqrtQ)
end