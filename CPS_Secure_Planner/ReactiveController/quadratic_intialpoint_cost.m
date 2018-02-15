function cost = quadratic_intialpoint_cost(x)
global block catXc terminal nx nu xO
xin = x-catXc;
xLast = xin(end-nx+1:end);
cost = x'*block*x + terminal*(norm(xLast(1:2)) + 1/norm(x(end-nx+1:end-nx+2)-xO) );
end

