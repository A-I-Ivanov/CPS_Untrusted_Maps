global deltaT block K nx nu catXc controlIndex xStart
T = 5;
deltaT = 0.1;
K = T/deltaT;
nx = 3;
nu = 2;

controlIndex = K*nx+1;
xC = [0,-10, -pi/2]';

catXc = repmat(xC, K,1);

R = [1 0 0; 0 1 0; 0 0 1];
Q = [1 0; 0 1];
terminal = 40*eye(3);

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

xStart = [0,0,pi/2]';

x0 = repmat(xStart, K,1);
x0 = vertcat(x0, zeros(nu*(K-1),1));

lb = zeros(nx*K+nu*(K-1),1);
ub = lb;
lb = lb-Inf;
ub = ub+Inf;
ub(controlIndex:2:end) = 2;
lb(controlIndex:2:end) = 0;
ub(controlIndex+1:2:end) = 1;
lb(controlIndex+1:2:end) = -1;


opts = optimset('Display','iter','Algorithm','interior-point', 'MaxIter', 100000, 'MaxFunEvals', 100000);
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);

[u,v] = pol2cart(result(3:3:controlIndex-1-nx),1);

quiver(result(1:3:controlIndex-1-nx), result(2:3:controlIndex-1-nx),u,v);




