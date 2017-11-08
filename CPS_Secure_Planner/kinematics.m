function xNew = kinematics(x,u)
global deltaT nx
xNew = zeros(nx,1);
xNew(1) = x(4)*cos(x(3))*deltaT + x(1);
xNew(2) = x(4)*sin(x(3))*deltaT + x(2);
xNew(3) = x(3) + x(5)*deltaT;
xNew(4) = x(4) +u(1)*deltaT;
xNew(5) = x(5) +u(2)*deltaT;
end