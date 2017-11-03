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

function xNew = kinematics(x,u)
global deltaT nx
xNew = zeros(nx,1);
xNew(1) = x(4)*cos(x(3))*deltaT + x(1);
xNew(2) = x(4)*sin(x(3))*deltaT + x(2);
xNew(3) = x(3) + x(5)*deltaT;
xNew(4) = x(4) +u(1)*deltaT;
xNew(5) = x(5) +u(2)*deltaT;
end


