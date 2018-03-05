%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cost, g] = minTimeCost(x)
global nx nu K xT

cost = x(1)*K;

g= zeros(length(x),1);
g(1) = K;
end


