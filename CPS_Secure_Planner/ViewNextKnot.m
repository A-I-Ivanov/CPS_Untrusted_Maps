%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
global deltaT  K nx nu controlIndex xStart xO rSafe rReact rSensor thetaSensor xT terminal
close all

T = 2;
deltaT = 0.1;
K = ceil(T/deltaT); %Number of time steps
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 5;
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .3;
thetaSensor = pi/2.1; %The angular range we care about (roughly in front of us)

controlIndex = K*nx+1; %I concatinate states and controls, this is for convenience
xT = [0,1.4, pi/2, .1, 0]'; %My referance point (the position I wnant to be close to)
xO = [0.01, rReact]';            %The location of my obstacle


R = [1 0 0 0 0;   %My state cost. Dont turn too much
    0 1 0 0 0;    %I don't care about the y so much
    0 0 1 0 0;     %How far off from the referance orientation is the robot?
    0 0 0 0 0;  %I don't care about velocities so they are 0
    0 0 0 0 0];

Q = [1 0; %%Control costs
    0 1];

terminal = 1000; %%This is a 'terminal cost'. 
               %How far away from the referance are you at the very end?

nonlcon = @knotViewConstraints;
fun = @minTimeCost;
A = [];
b = [];

xStart = [0,-.2,pi/2, .5, 0]'; %Where and how fast am I going?

x0 = xStart;
xlast = x0;

%The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]
uStart =  repmat([1 0]',(K-1),1);
for i =1:K-1
    xNext = kinematics(xlast,[1;0;]);
    xlast = xNext;
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end

x0 = vertcat(x0, uStart); %An initial guess of my controls too

Aeq = zeros(nx,length(x0));
Aeq(1:nx,1:nx) = eye(nx);
 
beq = xStart;


lb = zeros(nx*K+nu*(K-1),1);
ub = lb;
lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
ub = ub+Inf;
%ub(3:nx:controlIndex-1) = pi;
%lb(3:nx:controlIndex-1) = -pi;
lb(4:nx:controlIndex-1) = 0;
ub(5:nx:controlIndex-1) = pi/deltaT - .1;
lb(5:nx:controlIndex-1) = -pi/deltaT +.1;
ub(controlIndex:nu:end) = 2;  %%Bound the acceleration
lb(controlIndex:nu:end) = -2;    
ub(controlIndex+1:nu:end) = 10; %%Bound the turn rate
lb(controlIndex+1:nu:end) = -10;


%Some options for the optimizer
opts = optimset('Display','iter','Algorithm','sqp', 'MaxIter', 100000, 'MaxFunEvals', 100000,'TolX', 1e-14);
%opts.StepTolerance = 1e-10;
opts.ConstraintTolerance = 1e-5;

%What will my controller do for a range of object detections? I assume
%the obstale will be detected at the edge of rReact
sensLength = 2*cos(thetaSensor);

hold on


%%Solve the thing
tic
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);
toc

[u,v] = pol2cart(result(3:nx:controlIndex-1-nx),result(4:nx:controlIndex-1-nx));

%%Plot the x,y solutions
quiver(result(1:nx:controlIndex-1-nx), result(2:nx:controlIndex-1-nx),u,v,'b');
axis([-1 1 -1 1])

plot(xO(1),xO(2),'rx')
viscircles(xO',rSafe) %How far away did the robot need to be from the obstacle? Plot this


