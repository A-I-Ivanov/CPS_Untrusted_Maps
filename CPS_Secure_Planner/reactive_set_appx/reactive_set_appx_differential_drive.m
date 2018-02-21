%This code simulates a robot forward for 3 seconds using an MPC controller
clear all
close all
global deltaT  K nx nu controlIndex xStart xO rSafe rReact rSensor thetaSensor xT terminal
global num1 num2 block catXc 
num1 =0; num2=0;
figure

addpath('../ReactiveController');
deltaV = 0.1;
velDisc = 18;
elipseParams = zeros(velDisc, 4);

deltaT = 0.1;
K = 20; %Number of time steps
Tf = K*deltaT; %Final time
nx = 5; %The state number [x,y, orientation, linear speed, angular speed]
nu = 2; %The number of controls [acceleration, turn rate]
rSensor = 5;
rReact = .8; %The sensed range of an obstacle (when should our controller activate?)
rSafe = .5;
thetaSensor = pi/4; %The angular range we care about (roughly in front of us)


reactRange = 1; %When should the robot begin evasive maneuvers?
reactiveDisc = 10;

data = zeros(K*reactiveDisc, 2);


%Note when pointing straight down with one obstacle at [+-.0001, -rReact]'
%the solver fails. No known reason yet
controlIndex = K*nx+1; %I concatinate states and controls, this is for convenience
xT = [0,1.4, (pi/2), .2, 0]'; %My referance point (the position I wnant to be close to)
xO = [0, rReact]';            %The location of my obstacle


R = [0 0 0 0 0 0 0 0;   %Timing cost
     0 2 0 0 0 0 0 0;   %My state cost. Dont turn too much
     0 0 1 0 0 0 0 0;    %I don't care about the y so much
     0 0 0 0 0 0 0 0;     %How far off from the referance orientation is the robot?
     0 0 0 0 0 0 0 0;  %I don't care about velocities so they are 0
     0 0 0 0 0 0 0 0
     0 0 0 0 0 0 .1 0  %Control Cost 1
     0 0 0 0 0 0 0 .1]; %Control Cost 2

blkR = repmat({R},1,K);
block = blkdiag(blkR{:});

 
terminal = 10; %%This is a 'terminal cost'. 
               %How far away from the referance are you at the very end?

nonlcon = @kinematicWrapper;
fun = @quadratic_intialpoint_cost;
A = [];
b = [];


for k = 1:velDisc

    xStart = [0,-.2, pi/2, .1+k*deltaV, 0]'; %Where and how fast am I going?
    catXc = repmat([0; xStart(1:2); 0; 0; 0; 0; 0;], K, 1);
    uStart =  [-1 0]';

    x0 = [deltaT; xStart; uStart];

    xlast = x0;


    %The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

    for i =1:K-1
        stateOnly = diffDriveKinematics(xlast(2:1+nx),xlast(2+nx:end), xlast(1));
        xlast = [deltaT; stateOnly; uStart];
        x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
    end

    %x0 = vertcat(x0, uStart); %An initial guess of my controls too

    Aeq = zeros(nx+2,length(x0));
    Aeq(1:nx,2:nx+1) = eye(nx);
    Aeq(nx+1,:) = repmat([1 0 0 0 0 0 0 0], 1, K);
    Aeq(nx+2,end-nu-1) = 1;
    beq = [xStart; Tf; 0;];



    lb = zeros((nx+nu+1)*K,1);
    ub = lb;
    lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
    ub = ub+Inf;
    lb(1:nx+nu+1:end) = deltaT;
    ub(1:nx+nu+1:end) = deltaT;
    lb(5:nx+nu+1:end) = 0;
    %ub(5:nx:controlIndex-1) = pi/deltaT - .1;
    %lb(5:nx:controlIndex-1) = -pi/deltaT +.1;
    ub(nx+2:nx+nu+1:end) = 1;  %%Bound the acceleration
    lb(nx+2:nx+nu+1:end) = -1;    
    ub(nx+nu+1:nx+nu+1:end) = 4; %%Bound the turn rate
    lb(nx+nu+1:nx+nu+1:end) = -4;


    %Some options for the optimizer
    opts = optimset('Algorithm','sqp', 'MaxIter', 100000, 'MaxFunEvals', 100000,'TolX', 1e-16);
    %opts.StepTolerance = 1e-14;
    opts.ConstraintTolerance = 1e-6;

    %What will my controller do for a range of object detections? I assume
    %the obstale will be detected at the edge of rReact
    sensLength = reactRange*cos(thetaSensor);

    objectOffset = linspace(-sensLength, sensLength, reactiveDisc);%I assume
                                %the obstale will be detected at the edge of rS
    figure 
    hold on
    profile on
    for i=1:reactiveDisc
    xO = [objectOffset(i), reactRange]';
    %%Solve the thing
    
    result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);
    

    [u,v] = pol2cart(result(4:(nx+nu+1):end),result(5:(nx+nu+1):end));

    %%Plot the x,y solutions
    
    quiver(result(2:(nx+nu+1):end), result(3:(nx+nu+1):end),u,v,'b');
    axis([-1 1 -1 1])
    
    
    data((i-1)*(K)+1:i*(K), :) = [result(2:(nx+nu+1):end), result(3:(nx+nu+1):end)];

    %plot(xO(1),xO(2),'rx')
    %viscircles(xO',rSafe) %How far away did the robot need to be from the obstacle? Plot this

    %figure

    end
    
    elipseParams(k, 1:4) = boundingElipse(data(:,1), data(:,2));
    t = linspace(-pi,pi, 1000);
    elipseX = cos(t);
    elipseY = sin(t);

    elipseCoord = horzcat(elipseX',elipseY');

    Q = diag(sqrt(elipseParams(k, 3:4)));

    elipseCoord = elipseCoord*Q;

    plot(elipseCoord(:,1) + elipseParams(k,1), elipseCoord(:,2) + elipseParams(k,2),'g');

    
end


profile off
    profile viewer


