%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, constraints] = kinematicWrapper(x)
global nx nu K xO rSafe
    constraints = zeros(nx*(K-1),1);
    cin = zeros(K,1);
    deltaT = x(1);
    for i=1:K-1
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        constraints((i-1)*(nx)+1:i*nx) = x((nx+nu+1)*(i)+2:(nx+nu+1)*(i)+1+nx) - diffDriveKinematics(xNow, uNow, deltaT);
        cin(i) = -norm(xNow(1:2)-xO)^2+rSafe^2;
    end
    
    
end




