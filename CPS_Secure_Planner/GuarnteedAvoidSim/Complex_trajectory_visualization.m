%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
%close all
addpath('../BasicFunctions') ;
addpath('../ReactiveController');
global deltaT  K nx nu  xStart rSafe xT  polygons
global catXc turnBounds velBounds accelBounds unknownObst theFinal
global reactiveDT rReact rSensor  thetaSensor

numObst =1; %The number of obstacles 

%These variables are used to warmstart the algorithm or just look at 
%ouptuts by loading previously computed trajectories without solving 
testing =0;
warmstart = 0;

deltaT = .2; %Initial guess of deltaT
reactiveDT = .1; %DeltaT for the reactive controller 
K = 15; %Number of time steps or knot points
Tf = K*deltaT; %Final time
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 2; %Sensor range 
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .1; %Avoidance distance for the reactive controller
thetaSensor = pi/3; %The angular range of the sensor


%These are bounds for acceleration and turn rate as well as velocity
turnBounds = [-6 6];
velBounds = [0.01, 1.8];
accelBounds = [-1.6 1.6];

         
%Load known and unknown obstacles, initial condition, control guess, terminal condition
%To run the examples simply load the correct initial conditions using the
%files below.
% load('round_the_bend_test.mat');
%load('passage_way_large.mat');
%load('passage_way_small.mat');
% [a,b] = poly2cw(polygons{1}(1,:), polygons{1}(2,:));
% polygons{1} = [a;b];
thetaSensor = pi/3;
%%The Following is an example of how to set variables  to run your own tests
%%All of these variables are contained in the files above for the senarios
%%in our published work. 

xStart = [2.4 2.25 0 .7 0]'; %Where and how fast am I going?
xT = [3.1 1.625]'; %Terminal point 
theFinal = 0;
uStart =  [0 0]'; %guess of optimal controls
polygons{1} = [1 1 2.5 2.5;
                0 1 1 0];
polygons{2} = [1 1 2.5 2.5;
               1.5 2 2 1.5];
polygons{3} = [.5 .5 1.75 1.75;
               2.5 3 3 2.5];
polygons{4} = [3 3 3.5 3.5;
               1.75 3 3 1.75];
           
polygons{5} = [3 3 3.5 3.5;
               .5 1.5 1.5 .5];
           
           
           
TrajectoryAnimation(Traj,polygons, unknownObst, 1,1);
           