%%%%Written by Alexander I. Ivanov - 2017%%%%
function controls = reactiveControl(xStart, obstaclePoint, xSetpoint, newDT, rSteps, warmstart)

global  nx nu xO terminal turnBounds velBounds accelBounds reactiveSteps%top level globals
global  block catXc reactiveDT %this function's globals 

persistent Tf A b lb ub Aeq fun nonlcon options lastPath
if(isempty(Tf) || isempty(newDT)==0)
    
    if(isempty(newDT))
        reactiveDT = 0.1;
        reactiveSteps = 10; %Number of time steps
    else
        reactiveDT = newDT;
        reactiveSteps = rSteps;
    end
   
    Tf = reactiveSteps*reactiveDT; %Final time

    

    %The following cost assumes the robot is alligned with the y axis. 
    R = [ 2 0 0 0 0 0 0;   %My state cost. Dont turn too much
          0 2 0 0 0 0 0;    %I don't care about the y so much
          0 0 .001 0 0 0 0;     %How far off from the referance orientation is the robot?
          0 0 0 .001 0 0 0;  %I don't care about velocities so they are 0
          0 0 0 0 .001 0 0
          0 0 0 0 0 .1 0  %Control Cost 1
          0 0 0 0 0 0 .1]; %Control Cost 2
    blkR = [];
    blkR = horzcat(blkR,repmat({R},1,reactiveSteps-1));
    blkR = horzcat(blkR,{R(1:nx,1:nx)});
    block = blkdiag(blkR{:});

    totalVars = (nx+nu)*(reactiveSteps-1)+nx;
    terminal = 10; %%This is a 'terminal cost'. 
               %How far away from the referance are you at the very end?

    nonlcon = @reactiveControllerConst;
    fun = @quadratic_intialpoint_cost;
    A = [];
    b = [];

    catXc = [repmat([xSetpoint(1:2); 0; 0; 0; 0; 0;], reactiveSteps-1, 1);[xSetpoint(1:2); 0; 0; 0;];];

    Aeq = zeros(nx,totalVars);
    Aeq(1:nx,1:nx) = eye(nx);
    
    
    lb = zeros(totalVars,1);
    ub = lb;
    lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
    ub = ub+Inf;

    lb(4:nx+nu:end) = velBounds(1);
    ub(4:nx+nu:end) = velBounds(2);

    ub(nx+1:nx+nu:end) = accelBounds(2);  %%Bound the acceleration
    lb(nx+1:nx+nu:end) = accelBounds(1);    
    ub(nx+nu:nx+nu:end) = turnBounds(2); %%Bound the turn rate
    lb(nx+nu:nx+nu:end) = turnBounds(1);

    options = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', 200, 'MaxFunEvals', 20000,'TolX', 1e-16);
    options =optimoptions(options,'OptimalityTolerance', 1e-4);
    options =optimoptions(options,'ConstraintTolerance', 1e-5); 
    %options = optimoptions(options, 'UseParallel',true);
end

if(xStart(5)>0)
    uTurn = turnBounds(1);
else if(xStart(5)<0)
    uTurn = turnBounds(2);
    else
       uTurn =0; 
    end
end

uStart =  [accelBounds(1) uTurn]';


if(warmstart)
    x0 = [xStart ; lastPath(1+nx:end)];
else
     x0 = [xStart;];

    xlast = x0;

    stateOnly = xStart;
    %The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

    for i =1:reactiveSteps-1
        stateOnly = diffDriveKinematics(stateOnly,uStart, reactiveDT);
        xlast = [uStart; stateOnly];
        x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
    end
end



xO = obstaclePoint;            %The location of my obstacle

beq = [xStart;];


%%Solve the thing
tic
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
toc
lastPath = zeros(length(result)-nx-nu,1);
lastPath(1:end) = result(1+nx+nu:end);
controls = result(1+nx:nx+nu);    
  
  
end




