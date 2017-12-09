%This code simulates a robot forward for 3 seconds using an MPC controller
global deltaT catXc Sigmas nx nu TimeSteps terminal xT H Bdisc Adisc
close all
deltaT = 0.1;
T = 3;
nx = 2; %The state number [x,y, orientation, linear speed, angular speed]
nu = 1; %The number of controls [acceleration, turn rate]
nz = 3;
K = T/deltaT; %Number of time steps
TimeSteps =K;

A = [0 1;
     0 -.1;];
B = [0 .5;]';

odefun = @matrixExpTime;
[t,y] = ode45(odefun, [0, deltaT], [0;0;]);

Adisc = expm(A*deltaT);
Bdisc = y(end,:)';

H = [1 1 1; 0 0 0;]';

Sigma0 = eye(2);

SigmaV = [1 0 0; 0 2 0; 0 0 3;];

SigmaX = [.1 0; 0 .1;];

Sigmas = PreComputeSigmas(Sigma0, SigmaV, SigmaX, Adisc, H, K);

xStart = [-10;0]; %Where and how fast am I going?
xT = [-.1;0];

catXc = [zeros((K-1)*(nx+nu),1); xT;];

R = [1 0;   %My state cost. Dont turn too much
    0 .01;];    %Dont go too fast
    
Q = .01; %%Control costs


terminal = (K*eye(2)); %%This is a 'terminal cost'. 
                      %How far away from the referance are you at the very end?

block = R;
block = blkdiag(block, Q);

for i=1:K-2;
    block = blkdiag(block, R); %Make script R
    block = blkdiag(block, Q);
end

block = blkdiag(block, terminal);


nonlcon = [];
fun = @costCollision;
b = [];

DynamicConstr = [Adisc, Bdisc, -eye(nx)];
Aeq = [eye(2), zeros(size(DynamicConstr,1), size(block,1)-nx)];
Aeq = vertcat(Aeq,[DynamicConstr, zeros(size(DynamicConstr,1), size(block,1)-2*nx-nu)]);

for i =1:K-2
    head = zeros(size(DynamicConstr,1), (nx+nu)*i);
    tail = zeros(size(DynamicConstr,1), size(block,1)-size(DynamicConstr,2)-(nx+nu)*i);
    Aeq = vertcat(Aeq, [head, DynamicConstr, tail]);
end

Ain = zeros(1,size(Aeq,2));
Ain(end-1) = 1;
bin = 0;

beq = zeros(size(Aeq,1),1);
beq(1:2) = xStart;


%The optimizer needs a single vector. This is [x1, u1,....., xK,uK]
x0 = repmat([xT;0;], K-1,1);   %I neeed an initial guess of my trajecoty
x0 = vertcat(x0,xT);

lb = zeros(nx*K+nu*(K-1),1);
ub = lb;
lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
ub = ub+Inf;
ub(3:3:end) = 10;  %%Bound the acceleration
lb(3:3:end) = -10;    

%Some options for the optimizer :: 'Display','iter',
tic
opts = optimset('Algorithm','SQP', 'MaxIter', 100000, 'MaxFunEvals', 100000);
toc


%%Solve the thing
result = fmincon(fun,x0,Ain,bin,Aeq,beq,lb,ub,nonlcon, opts);


xhat = xStart;
xhatData = zeros(K,2);
xhatData(1,:) = xhat;
xData = zeros(K,2);
xData(1,:) = xStart';
uData = zeros(K-1,1);
a = zeros(3,1);
a(2) = -1;

for i =1:K-1
    uData(i) = result(1+nx);
    xData(i+1,:) = Adisc*xData(i,:)'+Bdisc*uData(i);
    z = H*xData(i+1,:)';
    
    if(i>0)
        z = z + a;
    end
    
    xhatData(i+1,:) = KalmanStep(xhatData(i,:)', uData(i), z, Sigmas{1}, SigmaX, SigmaV);
    
    
    beq = beq(3:end);
    beq(1:2) =  xhatData(i+1,:);
    Aeq = Aeq(1:size(Aeq,1)-nx, 1:end-(nx+nu));
    
    x0 = result(1+nx+nu:end);
    x0(1:2) = xhatData(i+1,:)';
    
    lb = lb(nx+nu+1:end);
    Sigmas = Sigmas(2:end);
    ub = ub(nx+nu+1:end);
    TimeSteps = TimeSteps-1;
    terminal = (TimeSteps*eye(2)); %%This is a 'terminal cost'. 
                      %How far away from the referance are you at the very end?
    
    Ain = zeros(1,size(Aeq,2));
    Ain(end-1) = 1;
    bin = 0;                  
                      
    result = fmincon(fun,x0,Ain,bin,Aeq,beq,lb,ub,nonlcon, opts);
end

%[u,v] = pol2cart(result(3:nx:controlIndex-1-nx),result(4:nx:controlIndex-1-nx));

%%Plot the x,y solutions
%quiver(result(1:nx:controlIndex-1-nx), result(2:nx:controlIndex-1-nx),u,v,'b');


Sigma0 = eye(2);

SigmaV = [1 0 0; 0 2 0; 0 0 3;];

SigmaX = [.1 0; 0 .1;];

Sigmas = PreComputeSigmas(Sigma0, SigmaV, SigmaX, Adisc, H, K);

hold on
plot(xData(:,1),linspace(0,T,K),'rx')
plot(xhatData(:,1),linspace(0,T,K),'bo')

hold off


title('Attack at t=0sec')
xlabel('Time sec')
ylabel('Position away from wall')
legend('Truth', 'Estimate', 'Location', 'northwest')


