function cost = quadratic_intialpoint_cost(x)
global block catXc
xin = x-catXc;
cost = xin'*block*xin;
end

