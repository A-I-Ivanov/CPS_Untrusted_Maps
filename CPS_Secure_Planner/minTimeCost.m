function cost = minTimeCost(x)
global K nx xT terminal
xX = x(1:nx:K*nx);
xY = x(2:nx:K*nx);
cost = sum((xX-xT(1)).^2 + (xY-xT(2)).^2) + terminal*((xX(end)-xT(1))^2 + (xY(end)-xT(2))^2);
end


