%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, constraints] = knotViewConstraints(x)
global nx nu K xStart xO rSafe rSensor rReact thetaSensor num2 num1
    a = [0.1978 0];
    Q = [0.0391 0; 0 0.0151];
    numConst = 1; 

    constraints = zeros(nx*(K-1),1);
    cin = zeros(K,1);
    deltaT = x(1);
    for i=1:K-2
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        constraints((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        
        x2C = x((nx+nu)*(i+1)+2:(nx+nu)*(i+1)+1+nx);
        
        cin(numConst*(i-1)+1) = -obstDistance(xNow(1:2))+rSafe;
        
        x2C(3) = wrapToPi(x2C(3)); %angle wrap
        xNow(3) = wrapToPi(xNow(3)); %angle wrap
        
        %The remaining constraints are safety visual contraints
        %They utilze eliptic transformations 
        if(x2C(4)<0)
            toaster=1;
        end
        [a,Q] = interpReactive(x2C(4));
        aRot = [cos(x2C(3)) -sin(x2C(3)); sin(x2C(3)) cos(x2C(3))]*a';
        xDiff = xNow(1:2)- x2C(1:2) - aRot;
        
        %%Has to be in sensor range, normalize to 1
        cin(numConst*(i-1)+2) = (norm(xDiff - aRot) - rSensor)/rSensor; 
        theta12 = -wrapToPi(atan2(xDiff(2), xDiff(1))) + pi;
        
        rotMat = [cos(theta12) -sin(theta12); sin(theta12) cos(theta12)]; 
        xTic = rotMat*xDiff;
        
        
        %Make sure we actually have intercepts
        if(wrapToPi(thetaSensor + xNow(3) + theta12)>(pi/2))
            m1 = 100;
            num1 = num1+1;
        else
            m1 = tan(thetaSensor + xNow(3) + theta12);
        end
        
        if(wrapToPi(-thetaSensor + xNow(3) + theta12)<-(pi/2))
            m2 = -100;
            num2 = num2+1;
        else
            m2 = tan(-thetaSensor + xNow(3) + theta12);
        end
        

        yTic1 = -[0 m1*xTic(1)]';
        yTic2 = -[0 m2*xTic(1)]';
        
        
        QTic = rotMat*Q*rotMat';
        
        sqrtQTic = chol(inv(QTic));
        
        x2Tic = sqrtQTic*xTic;
        y2Tic1 = sqrtQTic*yTic1;
        y2Tic2 = sqrtQTic*yTic2;
        
        A1 = (y2Tic1-x2Tic);
        A2 = (y2Tic2-x2Tic);
        b = -x2Tic;
        
        %Also these should be in linear inequality constriants
        %We normalize to one to get good constraint conditioning
        cin(numConst*(i-1)+3) = (-norm(A1*pinv(A1)*b + x2Tic) + 1)/(norm(x2Tic));
        cin(numConst*(i-1)+4) = (-norm(A2*pinv(A2)*b + x2Tic) + 1)/(norm(x2Tic));
        cin(numConst*(i-1)+5) = -y2Tic1(2);
        cin(numConst*(i-1)+6) = y2Tic2(2);
        
        
    end
    
        for i=K-1
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        constraints((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        cin(numConst*(i-1)+1) = -obstDistance(xNow(1:2))+rSafe;
        end
    
end



