%%%%Written by Alexander I. Ivanov - 2017%%%%
function cost = minTimeCost(x)
global nx nu K xT
terminal = 1000;
cost = x(1)*K;% + norm(x(end-nx-nu+1:end-nx-nu+2) - xT(1:2))*terminal ;
end


