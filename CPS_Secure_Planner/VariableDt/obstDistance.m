%%%%Written by Alexander I. Ivanov - 2017%%%%
%%This is a simple distnace function to point obstacles
function dist =  obstDistance(x)
global xO
numObst = size(xO);
 %Should do this with KDtree
 dist = Inf;
 for i=1:numObst(1)
     distNow = norm(x(1:2)'-xO(i,:));
    if( distNow < dist)
     dist = distNow;
    end
 end
 
end