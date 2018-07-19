function f = diffDriveDynamics(x,u)
f= zeros(5,1);

f(1) = x(4)*cos(x(3));
f(2) = x(4)*sin(x(3));
f(3) = x(5);
f(4) = u(1);
f(5) = u(2);

end