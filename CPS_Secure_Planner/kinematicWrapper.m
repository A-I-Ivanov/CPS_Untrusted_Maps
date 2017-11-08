function [cin, constraints] = kinematicWrapper(x)
global controlIndex nx nu K xStart xO rSafe
    constraints = zeros(nx*K,1);
    constraints(1:nx) = x(1:nx) - xStart;
    for i=1:K-1
        xNow = x(nx*(i-1)+1:nx*i);
        uNow = x(controlIndex+nu*(i-1):nu*i+controlIndex-1);
        constraints(nx*i+1:nx*(i+1)) = x(nx*i+1:nx*(i+1)) - kinematics(xNow, uNow);
        cin(i) = -norm(xNow(1:2)-xO)^2+rSafe^2;
    end
    
    
end




