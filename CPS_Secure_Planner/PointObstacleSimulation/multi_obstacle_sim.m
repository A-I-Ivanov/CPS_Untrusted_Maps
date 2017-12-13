%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
close all
global deltaT  K nx nu  xStart xO rSafe rReact rSensor thetaSensor xT 
global num1 num2  catXc
num1 =0; num2=0;

deltaT = 0.1;
K = 20; %Number of time steps
Tf = K*deltaT; %Final time
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 5;
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .3;
thetaSensor = pi/4; %The angular range we care about (roughly in front of us)



%Note when pointing straight down with one obstacle at [+-.0001, -rReact]'
%the solver fails. No known reason yet
xT = [0,2, (pi/2), .2, 0]'; %My referance point (the position I wnant to be close to)
xO = [0, rReact]';            %The location of my obstacle



nonlcon = @knotViewConstraints;
fun = @minTimeCost;
A = [];
b = [];

xStart = [0,-.3, pi/2, 1, 0]'; %Where and how fast am I going?
catXc = repmat([0; xStart(1:2); 0; 0; 0; 0; 0;], K, 1);
uStart =  [0 0]';

x0 = [deltaT; xStart; uStart];

xlast = [xStart; uStart];


%The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

for i =1:K-1
    stateOnly = diffDriveKinematics(xlast(1:nx),uStart, x0(1));
    xlast = [stateOnly; uStart];
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end


Aeq = zeros(2*nx,length(x0));
Aeq(1:nx,2:nx+1) = eye(nx);
Aeq(nx+1:2*nx,end-nu-nx+1:end-nu) = eye(nx);
beq = [xStart; xT;];



lb = zeros((nx+nu)*K +1,1);
ub = lb;
lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
ub = ub+Inf;
lb(1) = .001;
ub(1) = .2;
lb(5:nx+nu:end) = 0;
%ub(5:nx:controlIndex-1) = pi/deltaT - .1;
%lb(5:nx:controlIndex-1) = -pi/deltaT +.1;
ub(nx+2:nx+nu:end) = 1;  %%Bound the acceleration
lb(nx+2:nx+nu:end) = -1;    
ub(nx+nu+1:nx+nu:end) = 4; %%Bound the turn rate
lb(nx+nu+1:nx+nu:end) = -4;


%Some options for the optimizer
opts = optimset('Display','iter','Algorithm','sqp', 'MaxIter', 100000, 'MaxFunEvals', 100000, 'TolX', 1e-20);
%opts.StepTolerance = 1e-14;
opts.ConstraintTolerance = 1e-6;


figure 
hold on


%%Solve the thing
tic
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);
toc

[u,v] = pol2cart(result(4:(nx+nu):end),result(5:(nx+nu):end));

%%Plot the x,y solutions

quiver(result(2:(nx+nu):end), result(3:(nx+nu):end),u,v,'b');
axis([-1 1 -1 1])



plot(xO(1),xO(2),'rx')
viscircles(xO',rSafe) %How far away did the robot need to be from the obstacle? Plot this




