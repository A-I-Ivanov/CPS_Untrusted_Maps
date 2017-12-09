function cost = costWrapperSparse(x)
global block catXc
xin = x-catXc;
cost = xin'*block*xin;
end
