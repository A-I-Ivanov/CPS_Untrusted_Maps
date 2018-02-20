%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
close all
addpath('../BasicFunctions') ;
addpath('../ReactiveController');
global deltaT  K nx nu  xStart rSafe rReact rSensor thetaSensor xT  polygons
global num1 num2  catXc turnBounds velBounds accelBounds unknownObst
num1 =0; num2=0;

numObst =1;

testing =0;
warmstart = 0;

deltaT = 0.1;
reactiveDT = .1;
K = 15; %Number of time steps
Tf = K*deltaT; %Final time
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 2;
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .05;
thetaSensor = pi/3; %The angular range we care about (roughly in front of us)


turnBounds = [-4 4];
velBounds = [0.0, 1.8];
accelBounds = [-1.2 1.2];


%Note when pointing straight down with one obstacle at [+-.0001, -rReact]'
%the solver fails. No known reason yet
%xT = [1.2 , 1.5, (pi/2)]'; % .5, 0]'; %My referance point (the position I wnant to be close to)
            %The location of my obstacles
xT = [0 , 1.2, (pi/2)]'; % .5, 0]'; %My referance point (the position I wnant to be close to)
            
%Load test known polygons            
%load('large_box_known.mat');

polygons{1} = [-0.5,-.25,-.25;0,0,.5];
polygons{2} = [0.5,.25,.25;0,0,.5];


%Load unknown polygons
%load('small_box_unknown.mat');
                  

A = [];
b = [];

%xStart = [0.0, .1, 0, .5, 0]'; %Where and how fast am I going?
xStart = [0.0, -.5, pi/2, 1, 0]'; %Where and how fast am I going?
catXc = repmat([0; xStart(1:2); 0; 0; 0; 0; 0;], K, 1);
uStart =  [.1 0]';

x0 = [deltaT; xStart; uStart];



xlast = [xStart; uStart];


%The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

for i =1:K-2
    stateOnly = diffDriveKinematics(xlast(1:nx),uStart, x0(1));
    xlast = [stateOnly; uStart];
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end

stateOnly = diffDriveKinematics(xlast(1:nx),uStart, x0(1));
x0 = vertcat(x0, stateOnly);
%x0(8:end) = warmstart(8:end);
%x0(1) = warmstart(1);

Aeq = zeros(2*nx-2,length(x0));
Aeq(1:nx,2:nx+1) = eye(nx);
Aeq(nx+1:end,end-nx+1:end-nx+3) = eye(nx-2);
%Aeq(nx+1:nx+2,end-nu-nx+1:end-nu-nx+2) = eye(2);
beq = [xStart;xT]; % xT(1:2);];



lb = zeros((nx+nu)*(K-1)+1+nx,1);
ub = lb;
lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
ub = ub+Inf;
lb(1) = .001;
%ub(1) = .3;
lb(5:nx+nu:end-nx) = 0.1;
ub(5:nx+nu:end-nx) = 1.8;
%ub(5:nx:controlIndex-1) = pi/deltaT - .1;
%lb(5:nx:controlIndex-1) = -pi/deltaT +.1;
ub(nx+2:nx+nu:end-nx) = accelBounds(2);  %%Bound the acceleration
lb(nx+2:nx+nu:end-nx) = accelBounds(1);

lb(nx+nu+1:nx+nu:end-nx) = turnBounds(1);
ub(nx+nu+1:nx+nu:end-nx) = turnBounds(2); %%Bound the turn rate



%Some options for the optimizer
options = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', 200, 'MaxFunEvals', 10000, 'TolX', 1e-16);
options =optimoptions(options,'OptimalityTolerance', 1e-6);
options =optimoptions(options,'Display','iter');
options =optimoptions(options,'ConstraintTolerance', 1e-6); 
options =optimoptions(options,'FiniteDifferenceType', 'central'); 
options =optimoptions(options,'GradObj', 'on'); 
options =optimoptions(options,'SpecifyObjectiveGradient',true);
options =optimoptions(options,'CheckGradients',false);
options =optimoptions(options,'SpecifyConstraintGradient', true);
options = optimoptions(options, 'UseParallel',true);



%%Solve the thing
profile on

nonlcon = @basicDynamicsConstraints;%knotViewConstraints; %
fun = @minTimeCost;


tic

if(testing)
    load('testingTraj.mat')
    if(warmstart)
        x0 = plannedTraj;
        plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
    end
else
    plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
end
toc

plotObstacles_Traj(plannedTraj,polygons, unknownObst,xStart, xT, 0);

%[realizedTraj, setPoint] = simulateRobot(plannedTraj);

%plotObstacles_Traj(realizedTraj,polygons, unknownObst,xStart, xT, 1,0, setPoint);


%%Now do the same thing with my method
nonlcon = @knotViewConstraints; %basicDynamicsConstraints;%
fun = @minTimeCost;
tic
if(testing || warmstart)
    load('guaranteedTraj.mat')
    if(warmstart)
        x0 = plannedTraj;
        plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
    end
else
    plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
end
toc

plotObstacles_Traj(plannedTraj,polygons, unknownObst,xStart, xT, 0);

%[realizedTraj, setPoint] = simulateRobot(plannedTraj);

%plotObstacles_Traj(realizedTraj,polygons, unknownObst,xStart, xT, 1,0, setPoint);


profile viewer



