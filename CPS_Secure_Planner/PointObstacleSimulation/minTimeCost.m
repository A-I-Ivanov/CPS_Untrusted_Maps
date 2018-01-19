function [cost, g] = minTimeCost(x)
global nx nu K xT
terminal = 1000;
cost = x(1)*K;% + norm(x(end-nx-nu+1:end-nx-nu+2) - xT(1:2))*terminal ;

g= zeros(length(x),1);
g(1) = K;
end


