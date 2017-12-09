function cost = collisionProbs(x)
global nx nu K Sigmas 
cost = zeros(K,1);
rt2 = sqrt(2);
for i =1:K
    cost(i) = (1+erf(x(i)/(rt2*Sigmas{i}(1,1))))/2;
end

end