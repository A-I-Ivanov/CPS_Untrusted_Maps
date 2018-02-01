function dist =  obstDistance(x)
global xO rSafe
numObst = size(xO);
 %Should do this with KDtree
 dist = Inf;
 for i=1:numObst(1)
     distNow = norm(x(1:2)'-xO(i,:)) - rSafe;
    if( distNow < dist)
     dist = distNow;
    end
 end
 
end