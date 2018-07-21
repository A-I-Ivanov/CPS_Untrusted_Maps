function df = diffDriveDynamicsDeriv(x,u)
df= zeros(5,7);

df(1,3) = -x(4)*sin(x(3));
df(1,4) = cos(x(3));

df(2,3) = x(4)*cos(x(3));
df(2,4) = sin(x(3));

df(3,5) = 1;

df(4,6) = 1;

df(5,7) = 1;

end