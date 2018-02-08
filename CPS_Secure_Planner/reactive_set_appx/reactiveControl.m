function controls = reactiveControl(xStart, obstaclePoint, xSetpoint)

global  nx nu xO rSafe terminal turnBounds velBounds accelBounds %top level globals
global  block catXc %this function's globals

nx = 5; nu =2;
turnBounds = [-4 4];
velBounds = [0.1, 1.8];
accelBounds = [-1 1];


persistent deltaT K Tf A b lb ub Aeq fun nonlcon     
if(isempty(Tf))
    deltaT = 0.1;
    K = 10; %Number of time steps
    Tf = K*deltaT; %Final time





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
    
    blkR = horzcat(blkR,repmat({R},1,K-1));
    blkR = horzcat(blkR,{R(1:nx,1:nx)});
    block = blkdiag(blkR{:});

    totalVars = (nx+nu)*(K-1)+1+nx;
    terminal = 10; %%This is a 'terminal cost'. 
               %How far away from the referance are you at the very end?

    nonlcon = @kinematicWrapper;
    fun = @quadratic_intialpoint_cost;
    A = [];
    b = [];

    catXc = [deltaT; repmat([xSetpoint(1:2); 0; 0; 0; 0; 0;], K-1, 1);[xSetpoint(1:2); 0; 0; 0;];];

    Aeq = zeros(nx+2,totalVars);
    Aeq(1:nx,2:nx+1) = eye(nx);
    Aeq(nx+1,:) = [1, repmat(zeros(1,nx+nu), 1, K-1), zeros(1, nx)];
    Aeq(nx+2,end-nu-1) = 1;
    
    lb = zeros(totalVars,1);
    ub = lb;
    lb = lb-Inf; %What are the upper and lower bounds of my states and controls?
    ub = ub+Inf;
    lb(1:nx+nu+1:end) = deltaT;
    ub(1:nx+nu+1:end) = deltaT;
    lb(5:nx+nu+1:end) = velBounds(1);
    ub(5:nx+nu+1:end) = velBounds(2);

    ub(nx+2:nx+nu+1:end) = accelBounds(2);  %%Bound the acceleration
    lb(nx+2:nx+nu+1:end) = accelBounds(1);    
    ub(nx+nu+1:nx+nu+1:end) = turnBounds(2); %%Bound the turn rate
    lb(nx+nu+1:nx+nu+1:end) = turnBounds(1);


end




uStart =  [-1 0]';

x0 = [deltaT; xStart;];

xlast = x0;

stateOnly = xStart;
%The optimizer needs a single vector. This is [x1, x2,....., u1,u2.....]

for i =1:K-1
    stateOnly = diffDriveKinematics(stateOnly,uStart, deltaT);
    xlast = [uStart; stateOnly];
    x0 = vertcat(x0, xlast);   %I neeed an initial guess of my trajecoty
end



xO = obstaclePoint;            %The location of my obstacle

beq = [xStart; Tf; 0;];


%Some options for the optimizer
opts = optimset('Algorithm','sqp','Display', 'Iter', 'MaxIter', 100000, 'MaxFunEvals', 100000,'TolX', 1e-17);
opts.ConstraintTolerance = 1e-5;

%%Solve the thing

result = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon, opts);

   
controls = result(nx+1:nx+nu+1);    
  
  
end




