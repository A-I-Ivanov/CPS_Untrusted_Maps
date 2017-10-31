global deltaT block K nx nu catXc controlIndex xStart
T = 4;
deltaT = 0.05;
K = T/deltaT;
nx = 5;
nu = 2;
obst = [0,2];

controlIndex = K*nx+1;
xC = [0,0, -pi/2, 0, 0]';

catXc = repmat(xC, K,1);

R = [5 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0];
Q = [1 0; 0 1];
terminal = blkdiag(40*eye(3), zeros(2));

block = R;

for i=1:K-2;
    block = blkdiag(block, R);
end

block = blkdiag(block, terminal);

for i=1:K-1;
    block = blkdiag(block, Q);
end


nonlcon = @kinematicWrapper;
fun = @costWrapper;
A = [];
b = [];
Aeq = [];
beq = [];

xStart = [0,0,pi/2, 1, 0]';

x0 = repmat(xStart, K,1);
x0 = vertcat(x0, zeros(nu*(K-1),1));

lb = zeros(nx*K+nu*(K-1),1);
ub = lb;
lb = lb-Inf;
ub = ub+Inf;
lb(4:nx:end) = 0;
ub(controlIndex:2:end) = 2;
lb(controlIndex:2:end) = -2;
ub(controlIndex+1:2:end) = 1;
lb(controlIndex+1:2:end) = -1;


opts = optimset('Display','iter','Algorithm','interior-point', 'MaxIter', 100000, 'MaxFunEvals', 100000);
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);

[u,v] = pol2cart(result(3:nx:controlIndex-1-nx),result(4:nx:controlIndex-1-nx));

quiver(result(1:nx:controlIndex-1-nx), result(2:nx:controlIndex-1-nx),u,v);




