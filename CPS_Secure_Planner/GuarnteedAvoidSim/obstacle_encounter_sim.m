%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
close all
addpath('../BasicFunctions') ;
addpath('../ReactiveController');
global deltaT  K nx nu  xStart rSafe xT  polygons
global catXc turnBounds velBounds accelBounds unknownObst
global reactiveDT rReact rSensor  thetaSensor

numObst =1; %The number of obstacles 

%These variables are used to warmstart the algorithm or just look at 
%ouptuts by loading previously computed trajectories without solving 
testing =0;
warmstart = 0;

deltaT = 0.2; %Initial guess of deltaT
reactiveDT = .1; %DeltaT for the reactive controller 
K = 20; %Number of time steps or knot points
Tf = K*deltaT; %Final time
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 2; %Sensor range 
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .1; %Avoidance distance for the reactive controller
thetaSensor = pi/3; %The angular range of the sensor


%These are bounds for acceleration and turn rate as well as velocity
turnBounds = [-4 4];
velBounds = [0.01, 1.8];
accelBounds = [-1.2 1.2];

         
%Load known and unknown obstacles, initial condition, control guess, terminal condition
%To run the examples simply load the correct initial conditions using the
%files below.
load('round_the_bend_test.mat');
%load('passage_way_large.mat');
%load('passage_way_small.mat');
[a,b] = poly2cw(polygons{1}(1,:), polygons{1}(2,:));
polygons{1} = [a;b];
thetaSensor = pi/4;
%%The Following is an example of how to set variables  to run your own tests
%%All of these variables are contained in the files above for the senarios
%%in our published work. 

%xStart = [0.0, -.5, pi/2, 1, 0]'; %Where and how fast am I going?
%xT = [0 , 1.2, (pi/2)]'; % .5, 0]'; %Terminal point 
uStart =  [.2 .4]'; %guess of optimal controls
%polygons{1} = [-0.5,-.25,-.25;0,0,.5];
%polygons{2} = [0.5,.25,.25;0,0,.5];


%Load unknown polygons
%load('small_box_unknown.mat');
                  
%Linear inequality constraint matricies 
A = [];
b = [];

%The following matrix is used for the reatice controller
catXc = repmat([0; xStart(1:2); 0; 0; 0; 0; 0;], K, 1);

x0 = [deltaT; xStart; uStart];



xlast = [xStart; uStart];


%The optimizer needs a single vector. This is [deltaT, x1, u1, x2, u2,..., x_K]
for i =1:K-2
    stateOnly = diffDriveKinematics(xlast(1:nx),uStart, x0(1));
    xlast = [stateOnly; uStart];
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end

stateOnly = diffDriveKinematics(xlast(1:nx),uStart, x0(1));
x0 = vertcat(x0, stateOnly);

%The equality constraints
Aeq = zeros(2*nx-2,length(x0));
Aeq(1:nx,2:nx+1) = eye(nx);
Aeq(nx+1:end,end-nx+1:end-nx+3) = eye(nx-2);
beq = [xStart;xT]; %Make sure the start and end conditions are satisfied



lb = zeros((nx+nu)*(K-1)+1+nx,1);
ub = lb;
lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
ub = ub+Inf;
lb(1) = .001; %Lower bound deltaT
lb(5:nx+nu:end-nx) = velBounds(1); %Bound velocity 
ub(5:nx+nu:end-nx) = velBounds(2);

ub(nx+2:nx+nu:end-nx) = accelBounds(2);  %%Bound the acceleration
lb(nx+2:nx+nu:end-nx) = accelBounds(1);

lb(nx+nu+1:nx+nu:end-nx) = turnBounds(1);
ub(nx+nu+1:nx+nu:end-nx) = turnBounds(2); %%Bound the turn rate



%Some options for the optimizer
options = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', 300, 'MaxFunEvals', 10000, 'TolX', 1e-12);
options =optimoptions(options,'OptimalityTolerance', 1e-3);
options =optimoptions(options,'Display','iter');
options =optimoptions(options,'ConstraintTolerance', 1e-4); 
options =optimoptions(options,'FiniteDifferenceType', 'central'); 
options =optimoptions(options,'GradObj', 'on'); 
options =optimoptions(options,'SpecifyObjectiveGradient',true);
options =optimoptions(options,'CheckGradients',false);
options =optimoptions(options,'SpecifyConstraintGradient', true);
options = optimoptions(options, 'UseParallel',true);



%%Solve the thing
% profile on
% safetyDistance = .025;
% nonlcon = @(x) basicDynamicsConstraints(x, safetyDistance);%Constriants for base-line planner
% fun = @minTimeCost;
% 
% 
% tic
% 
% if(testing)
%     %Uncomment this line for testing and visualization. 
%     %load('testingTraj.mat') %Save your base-line trajectory to this file if you want
%     %to simply visualize the path
%     if(warmstart)
%         x0 = plannedTraj;
%         plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
%     end
% else
%     plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
% end
% toc
% 
% plotObstacles_Traj(plannedTraj,polygons, unknownObst,xStart, xT, 0);
% 
% [realizedTraj, setPoint] = simulateRobot(plannedTraj);
% 
% plotObstacles_Traj(realizedTraj,polygons, unknownObst,xStart, xT, 1,0, setPoint);


%%Now do the same thing with secure method
nonlcon = @knotViewConstraints; %Visibility constraints, dynamics, and reactive set constraints
fun = @minTimeCost;
warmstart=1;
tic
if(testing || warmstart)
    safetyDistance = .1;
    nonlcon = @(x) basicDynamicsConstraints(x, safetyDistance);%Constriants for base-line planner
    plannedTraj = fmincon(fun,plannedTraj,A,b,Aeq,beq,lb,ub,nonlcon, options);
    fun = @minTimeCost;
    nonlcon = @knotViewConstraints;
    %Uncomment this line for testing and visualization. 
    %load('guaranteedTraj.mat') %Save your optimal trajectory to this file if you want
    %to simply visualize the path
    if(warmstart) %Note: If the warmstart path is highly infeasible, this will fail
        plannedTraj = fmincon(fun,plannedTraj,A,b,Aeq,beq,lb,ub,nonlcon, options);
    end
else
    plannedTraj = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
end
toc

%To visualize all reactive set elipses uncomment the following line.
%plotObstacles_Traj(plannedTraj,polygons, unknownObst,xStart, xT, 0,1);
plotObstacles_Traj(plannedTraj,polygons, unknownObst,xStart, xT, 0);

[realizedTraj, setPoint] = simulateRobot(plannedTraj);

%To visualize all reactive set elipses uncomment the following line.
%plotObstacles_Traj(realizedTraj,polygons, unknownObst,xStart, xT, 1,1);
plotObstacles_Traj(realizedTraj,polygons, unknownObst,xStart, xT, 1,0, setPoint);


profile viewer



