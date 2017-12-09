function cost = costCollision(x)
global nx nu Timesteps Sigmas terminal xT
cost =0;
pos = x(1:(nx+nu):end);
rt2 = sqrt(2);
for i =1:Timesteps
    costk = log(-erf(pos(i)/(rt2*Sigmas{i}(1,1)))+1);
    if costk<-500
        costk = -500;
    end
    cost = cost - costk;
end

cost = cost+ (x(end-nx+1:end)-xT)'*terminal*(x(end-nx+1:end)-xT);
end