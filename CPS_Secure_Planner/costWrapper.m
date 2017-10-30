function cost = costWrapper(x)
global K nx
cost = QuadCost(x(1:K*nx), x(K*nx+1:end));
end

function cost = QuadCost(x, u)
global block catXc
xin = x-catXc;
cost = [xin;u]'*block*[xin;u];
end

