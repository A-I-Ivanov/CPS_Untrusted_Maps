%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cost, g] = minTimeCost(x)
global nx nu K xT theFinal
scale = 100;
theCost = 10;
controls1 = x(7:9:end)/scale;
controls2 = x(8:9:end)/scale;
controls3 = x(9:9:end)/scale;
controls4 = x(10:9:end)/scale;

cost = x(1)*x(1)*K/2 + (sum(controls1.^2)+sum(controls2.^2)+sum(controls3.^2)+sum(controls4.^2))/2;
cost = cost + theCost*(x(end-nu-3)-theFinal)*(x(end-nu-3)-theFinal)/2;
g= zeros(length(x),1);
g(1) = K*x(1);
g(7:9:end) = controls1/scale;
g(8:9:end) = controls2/scale;
g(9:9:end) = controls3/scale;
g(10:9:end) = controls4/scale;
g(end-nu-3) = theCost*(x(end-nu-3)-theFinal);
end


