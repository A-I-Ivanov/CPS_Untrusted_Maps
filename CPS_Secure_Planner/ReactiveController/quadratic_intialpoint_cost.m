%%%%Written by Alexander I. Ivanov - 2017%%%%
function cost = quadratic_intialpoint_cost(x)
global block catXc terminal nx nu xO
xin = x-catXc;
xLast = xin(end-nx+1:end);
poseX = x(1:nx+nu:end) - xO(1);
poseY = x(2:nx+nu:end) - xO(2);
cost = x'*block*x + terminal*(norm(xLast(1:2))) + 1/norm(x(end-nx+1:end-nx+2)-xO) ;
end

