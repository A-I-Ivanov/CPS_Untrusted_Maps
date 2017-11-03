global vectorized 

vectorized = [data(1:nx:controlIndex-1-nx,1), data(2:nx:controlIndex-1-nx,1)];
for i = 2:2*nplot
vectorized = vertcat(vectorized,[data(1:nx:controlIndex-1-nx,i), data(2:nx:controlIndex-1-nx,i)]);
end
hold off


fun = @(x) x(3)+x(4);

A =[];
b = [];
lb = zeros(4,1);
ub = lb;
lb(1:2) = -Inf;
ub = ub+Inf;
x0 = [0,0,.1,.1];
nonlcon = @elliptic_cost;
Aeq = [];
beq = [];


elipseParams = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);

t = linspace(-pi,pi, 1000);
elipseX = cos(t);
elipseY = sin(t);

elipseCoord = horzcat(elipseX',elipseY');

Q = diag(sqrt(elipseParams(3:4)));

elipseCoord = elipseCoord*Q;

hold on 
plot(elipseCoord(:,1) + elipseParams(1), elipseCoord(:,2) + elipseParams(2),'r');

