%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, constraints] = reactiveControllerConst(x)
global nx nu reactiveSteps xO rSafe deltaT
    constraints = zeros(nx*(reactiveSteps-1),1);
    cin = zeros(reactiveSteps,1);
    for i=1:reactiveSteps-1
        xNow = x((nx+nu)*(i-1)+1:(nx+nu)*(i-1)+nx);
        uNow = x((nx+nu)*(i-1)+nx+1:(nx+nu)*(i-1)+nu+nx);
        constraints((i-1)*(nx)+1:i*nx) = x((nx+nu)*(i)+1:(nx+nu)*(i)+nx) - diffDriveKinematics(xNow, uNow, deltaT);
        cin(i) = -norm(xNow(1:2)-xO)+rSafe;
    end
    
    cin(reactiveSteps) = -norm(x(end-nx+1:end-nx+2)-xO)+rSafe;
end




