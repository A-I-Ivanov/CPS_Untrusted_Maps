function cost = quadratic_intialpoint_cost(x)
global block catXc terminal nx nu xO
xin = x-catXc;
xLast = xin(end-nx-nu:end);
cost = xin'*block*xin + terminal*(norm(xLast(2:3)) - norm(x(end-nx-nu+1:end-nx-nu+2)-xO) );
end

