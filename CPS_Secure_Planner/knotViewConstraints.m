function [cin, constraints] = knotViewConstraints(x)
global controlIndex nx nu K xStart xO rSafe rSensor thetaSensor
    a = [0.1978 0];
    Q = [0.0391 0; 0 0.0151];
    
    constraints = zeros(nx*K,1);
    constraints(1:nx) = x(1:nx) - xStart;
    for i=1:K-1
        xNow = x(nx*(i-1)+1:nx*i);
        xNext = x(nx*i+1:nx*(i+1));
        uNow = x(controlIndex+nu*(i-1):nu*i+controlIndex-1);
        constraints(nx*i+1:nx*(i+1)) =  xNext - kinematics(xNow, uNow);
        cin(i) = -norm(xNow(1:2)-xO)^2+rSafe^2;
        
        %The remaining constraints are safety visual contraints
        %They utilze eliptic transformations 
        aRot = [cos(xNow(3)) -sin(xNow(3)); sin(xNow(3)) cos(xNow(3))]*a';
        xDiff = xNow(1:2)- xNext(1:2) - aRot;
        cin(i+K-1) = norm(xDiff - aRot) - rSensor-10; %%Has to be in sensor range
        theta12 = atan2(xDiff(2), xDiff(1));
        
        rotMat = [cos(-theta12) -sin(-theta12); sin(-theta12) cos(-theta12)]; 
        xTic = rotMat*xDiff;
        
        m1 = tan(thetaSensor - xNow(3) - theta12);
        m2 = tan(-thetaSensor - xNow(3) - theta12);
        
        yTic1 = [0 m1*xTic(1)]';
        yTic2 = [0 m2*xTic(1)]';
        
        
        QTic = rotMat*Q*rotMat';
        
        try
        sqrtQTic = chol(inv(QTic));
        catch
            QTic = Q;
            
        end
        
        x2Tic = sqrtQTic*xTic;
        y2Tic1 = sqrtQTic*yTic1;
        y2Tic2 = sqrtQTic*yTic2;
        
        A1 = (y2Tic1-x2Tic);
        A2 = (y2Tic2-x2Tic);
        b = -x2Tic;
        
        
        if norm(inv(A1'*A1)*A1'*b)>=1 || norm(inv(A2'*A2)*A2'*b)>=1
            gay=1;
        end
        cin(i+2*(K-1)) = -norm(A1*inv(A1'*A1)*A1'*b) + 1;
        cin(i+3*(K-1)) = -norm(A2*inv(A2'*A2)*A2'*b) + 1;
        
        
        
    end
    
    
end

function xNew = kinematics(x,u)
global deltaT nx
xNew = zeros(nx,1);
xNew(1) = x(4)*cos(x(3))*deltaT + x(1);
xNew(2) = x(4)*sin(x(3))*deltaT + x(2);
xNew(3) = x(3) + x(5)*deltaT;
xNew(4) = x(4) +u(1)*deltaT;
xNew(5) = x(5) +u(2)*deltaT;
end


