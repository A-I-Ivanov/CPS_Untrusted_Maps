function elipseParams = boundingElipse(xData, yData)
%This function fits the minimu bounding 2D elispe to a set of trajectory
%data
global vectorized 
    
    vectorized = [xData, yData];

    fun = @(x) x(3)+x(4);

    A =[];
    b = [];
    lb = zeros(4,1);
    ub = lb;
    lb(1:2) = -Inf;
    ub = ub+Inf;
    x0 = [0,0,.1,.1];
    nonlcon = @elliptic_constr;
    Aeq = [];
    beq = [];


    elipseParams = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
  
    
    
end
