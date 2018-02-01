function [cin, eqconst, DCIn, DCeq] = knotViewConstraints(x)
    finDiffDelta = 1e-8;
global nx nu K xStart xO rSafe rSensor rReact thetaSensor num2 num1

    E = eye(nx);
    numConst = 6; 
    numEqConst = nx*(K-1);
    deltaEye = 2;
    eqconst = zeros(numEqConst,1);
    
    DCeq = zeros(numEqConst, length(x)); 
    DCIn = zeros(numConst*(K-deltaEye), length(x)); %+deltaEye
    

    
    cin = zeros(numConst*(K-deltaEye),1);
    deltaT = x(1);
    
    %Define state-based constraints and visibility constraints
    for i=1:K-deltaEye-1
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        eqconst((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        %%Populate sparse jacobian matrix for differential drive
        DCeq((i-1)*(nx)+1:i*nx,:)=popSparseJac(i, DCeq((i-1)*(nx)+1:i*nx,:), xNow, uNow, deltaT);
        
        
        x2C = x((nx+nu)*(i+1)+2:(nx+nu)*(i+1)+1+nx);
        
        cin(numConst*(i-1)+1:numConst*(i)) = calcInEq(xNow, x2C);
        DCIn(numConst*(i-1)+1:numConst*(i),:) = popSparseInJac(i,DCIn(numConst*(i-1)+1:numConst*(i),:), cin(numConst*(i-1)+1:numConst*(i)), xNow, x2C);
        
        
    end
    
    %The following are just state-based constraints
    l = 0; %this makes other indexing easier
    for i=K-deltaEye:K-1
        l = l+1;
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        eqconst((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        %%Populate sparse jacobian matrix for differential drive
        DCeq((i-1)*(nx)+1:i*nx,:)=popSparseJac(i, DCeq((i-1)*(nx)+1:i*nx,:), xNow, uNow, deltaT);
        
        %cin(numConst*(K-deltaEye)+l) = -obstDistance(xNow);
        
        %Take care of terminal inequality constraints
        for k=1:nx
            vDelt = finDiffDelta*E(:,k);
            deltDist = -obstDistance(xNow +vDelt); %Check distance to obstacles
            %DCIn(140-deltaEye+l,174-((deltaEye-l+1)*(nx+nu))+nu+k) = (deltDist - cin(numConst*(K-deltaEye)+l))/finDiffDelta;
        end
        
        
    end
    
    DCIn = DCIn';
    DCeq = DCeq';
    
    return
    
    
    
    %%%%%%%Begin Helper Functions%%%%%%%%%%%
    
    %This function populates the sparse inequality jacobian using finite
    %differences
    function rows = popSparseInJac(i,DinRows, cin, xNow, x2C)
        rows = DinRows;
        for j=1:nx
            vectDelt = finDiffDelta*E(:,j);
            cDelt = calcInEq(xNow + vectDelt,x2C)';
            rows(:,1+j+(i-1)*(nx+nu))= (cDelt-cin)/finDiffDelta;
            cDelt = calcInEq(xNow,x2C  + vectDelt)';
            rows(:,1+j+(i-1+deltaEye)*(nx+nu))= (cDelt-cin)/finDiffDelta;
        end

    end
    
    %This function populates rows of the sparse equality jacobian analytically
    function rows = popSparseJac(i, DeqRows, xNow, uNow, deltaT)
        C = cos(xNow(3));
        S = sin(xNow(3));
        V = xNow(4);
        rows = DeqRows;
        
        rows(1,1) = -V*C;        
        rows(1,(i-1)*(nx+nu)+4) = V*S*deltaT;
        rows(1,(i-1)*(nx+nu)+5) = -C*deltaT;
        rows(1,(i-1)*(nx+nu)+9) = 1;
        
        rows(2,1) = -V*S;        
        rows(2,(i-1)*(nx+nu)+4) = -V*C*deltaT;
        rows(2,(i-1)*(nx+nu)+5) = -S*deltaT;
        rows(2,(i-1)*(nx+nu)+10) = 1;
        
        rows(3,1) = -xNow(5);        
        rows(3,(i-1)*(nx+nu)+6) = -deltaT;
        rows(3,(i-1)*(nx+nu)+11) = 1;
        
        rows(4,1) = -uNow(1);        
        rows(4,(i-1)*(nx+nu)+7) = -deltaT;
        rows(4,(i-1)*(nx+nu)+12) = 1;
        
        rows(5,1) = -uNow(2);        
        rows(5,(i-1)*(nx+nu)+8) = -deltaT;
        rows(5,(i-1)*(nx+nu)+13) = 1;
    end


    function ineqs = calcInEq(xNow, x2C)
        ineqs(1) = -obstDistance(xNow); %Check distance to obstacles
        
        x2C(3) = wrapToPi(x2C(3)); %angle wrap
        xNow(3) = wrapToPi(xNow(3)); %angle wrap
        
        %The remaining eqconst are safety visual contraints
        %They utilze eliptic transformations 
        [a,Q] = interpReactive(x2C(4));
        aRot = [cos(x2C(3)) -sin(x2C(3)); sin(x2C(3)) cos(x2C(3))]*a';
        xDiff = xNow(1:2)- x2C(1:2) - aRot;
        
        %%Has to be in sensor range, normalize to 1
        ineqs(2) = (norm(xDiff - aRot) - rSensor)/rSensor; 
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
        ineqs(3) = (-norm(A1*pinv(A1)*b + x2Tic) + 1)/(norm(x2Tic));
        ineqs(4) = (-norm(A2*pinv(A2)*b + x2Tic) + 1)/(norm(x2Tic));
        ineqs(5) = -y2Tic1(2);
        ineqs(6) = y2Tic2(2);
    end
    
end



