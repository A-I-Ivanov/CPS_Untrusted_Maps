function [cin, constraints] = knotViewConstraints(x)
global controlIndex nx nu K xStart xO rSafe rSensor rReact thetaSensor
    a = [0.1978 0];
    Q = [0.0391 0; 0 0.0151];
    numConst = 6; 
    constraints = zeros(nx*K,1);
    for i=1:K-1
        xNow = x(nx*(i-1)+1:nx*i);
        
        xNext = x(nx*i+1:nx*(i+1));
        
        uNow = x(controlIndex+nu*(i-1):nu*i+controlIndex-1);
        constraints(nx*i+1:nx*(i+1)) =  xNext - kinematics(xNow, uNow);
        cin(numConst*(i-1)+1) = -norm(xNow(1:2)-xO)+rSafe;
        
        xNext(3) = mod(xNext(3), 2*pi); %angle wrap
        xNow(3) = mod(xNow(3), 2*pi); %angle wrap
        
        %The remaining constraints are safety visual contraints
        %They utilze eliptic transformations 
        aRot = [cos(xNext(3)) -sin(xNext(3)); sin(xNext(3)) cos(xNext(3))]*a';
        xDiff = +xNow(1:2)- xNext(1:2) - aRot;
        
        %%Has to be in sensor range, normalize to 1
        cin(numConst*(i-1)+2) = (norm(xDiff - aRot) - rSensor)/rSensor; 
        theta12 = -mod(atan2(xDiff(2), xDiff(1)), 2*pi) + pi;
        
        rotMat = [cos(theta12) -sin(theta12); sin(theta12) cos(theta12)]; 
        xTic = rotMat*xDiff;
        
        if(xTic(1)>0)
            toaster=1;
        end
        
        %Make sure we actually have intercepts
        if(mod(thetaSensor + xNow(3) + theta12,2*pi)>(pi/2))
            m1 = 100;
        else
            m1 = tan(thetaSensor + xNow(3) + theta12);
        end
        
        if(mod(-thetaSensor + xNow(3) + theta12, 2*pi)>(3*pi/4))
            m2 = -100;
        else
            m2 = tan(-thetaSensor + xNow(3) + theta12);
        end
        
        if(m1<0)
            toaster=1;
        end
        
        if(m2>0)
            toaster=1;
        end

        
        yTic1 = -[0 m1*xTic(1)]';
        yTic2 = -[0 m2*xTic(1)]';
        
        
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
        %Poor use of sparsity here. Should re-factor
        %Also these should be in linear inequality constriants
        %We normalize to one to get good constraint conditioning
        cin(numConst*(i-1)+3) = (-norm(A1*inv(A1'*A1)*A1'*b + x2Tic) + 1)/(norm(x2Tic));
        cin(numConst*(i-1)+4) = (-norm(A2*inv(A2'*A2)*A2'*b + x2Tic) + 1)/(norm(x2Tic));
        cin(numConst*(i-1)+5) = -y2Tic1(2);
        cin(numConst*(i-1)+6) = y2Tic2(2);
        
        
        
    end
    
    if(max(cin)> 1e-6)
     a = max(cin)
     
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


