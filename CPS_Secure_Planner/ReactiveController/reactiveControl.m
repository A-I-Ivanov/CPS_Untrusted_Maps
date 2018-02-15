function controls = reactiveControl(xStart, obstaclePoint, xSetpoint, newDT, rSteps)

global  nx nu xO rSafe terminal turnBounds velBounds accelBounds reactiveSteps%top level globals
global  block catXc %this function's globals

rSafe = .3;
nx = 5; nu =2;
turnBounds = [-4 4];
velBounds = [0.0, 1.8];
accelBounds = [-1 1];


persistent deltaT Tf A b lb ub Aeq fun nonlcon options   
if(isempty(Tf) || isempty(newDT)==0)
    
    if(isempty(newDT))
        deltaT = 0.1;
        reactiveSteps = 10; %Number of time steps
    else
        deltaT = newDT;
        reactiveSteps = rSteps;
    end
   
    Tf = reactiveSteps*deltaT; %Final time

    %Note when pointing straight down with one obstacle at [+-.0001, -rReact]'
    %the solver fails. No known reason yet


    R = [ 2 0 0 0 0 0 0;   %My state cost. Dont turn too much
          0 1 0 0 0 0 0;    %I don't care about the y so much
          0 0 0 0 0 0 0;     %How far off from the referance orientation is the robot?
          0 0 0 0 0 0 0;  %I don't care about velocities so they are 0
          0 0 0 0 0 0 0
          0 0 0 0 0 .1 0  %Control Cost 1
          0 0 0 0 0 0 .1]; %Control Cost 2
    
    blkR{1} = 0;
    
    blkR = horzcat(blkR,repmat({R},1,reactiveSteps-1));
    blkR = horzcat(blkR,{R(1:nx,1:nx)});
    block = blkdiag(blkR{:});

    totalVars = (nx+nu)*(reactiveSteps-1)+1+nx;
    terminal = 10; %%This is a 'terminal cost'. 
               %How far away from the referance are you at the very end?

    nonlcon = @reactiveControllerConst;
    fun = @quadratic_intialpoint_cost;
    A = [];
    b = [];

    catXc = [deltaT; repmat([xSetpoint(1:2); 0; 0; 0; 0; 0;], reactiveSteps-1, 1);[xSetpoint(1:2); 0; 0; 0;];];

    Aeq = zeros(nx+1,totalVars);
    Aeq(1:nx,2:nx+1) = eye(nx);
    Aeq(nx+1,:) = [1, repmat(zeros(1,nx+nu), 1, reactiveSteps-1), zeros(1, nx)];
    
    lb = zeros(totalVars,1);
    ub = lb;
    lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
    ub = ub+Inf;
    lb(1) = deltaT;
    ub(1) = deltaT;
    lb(5:nx+nu:end) = velBounds(1);
    ub(5:nx+nu:end) = velBounds(2);

    ub(nx+2:nx+nu:end) = accelBounds(2);  %%Bound the acceleration
    lb(nx+2:nx+nu:end) = accelBounds(1);    
    ub(nx+nu+1:nx+nu:end) = turnBounds(2); %%Bound the turn rate
    lb(nx+nu+1:nx+nu:end) = turnBounds(1);

    options = optimoptions('fmincon','Algorithm','sqp', 'MaxIter', 200, 'MaxFunEvals', 20000,'TolX', 1e-16);
    options =optimoptions(options,'OptimalityTolerance', 1e-4);
    options =optimoptions(options,'ConstraintTolerance', 1e-5); 
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

x0 = [deltaT; xStart;];

xlast = x0;

stateOnly = xStart;
%The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

for i =1:reactiveSteps-1
    stateOnly = diffDriveKinematics(stateOnly,uStart, deltaT);
    xlast = [uStart; stateOnly];
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end



xO = obstaclePoint;            %The location of my obstacle

beq = [xStart; deltaT;];


%%Solve the thing
tic
result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
toc
   
controls = result(2+nx:nx+nu+1);    
  
  
end




